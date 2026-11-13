#include <chargefw/adapters/native/json_output.h>

#include "adapters/component.h"
#include "adapters/record_assignments.h"

#include <chargefw/version.h>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <optional>
#include <ostream>
#include <print>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace chargefw::adapters::native::json_output {
namespace {

using Json = nlohmann::json;

constexpr auto metric_scale = 1000.0;

[[nodiscard]] auto diagnostics_json(const std::span<const ResultDiagnostic> diagnostics) -> Json {
    const auto severity_name = [](const DiagnosticSeverity severity) -> const char* {
        switch (severity) {
        case DiagnosticSeverity::info:
            return "info";
        case DiagnosticSeverity::warning:
            return "warning";
        case DiagnosticSeverity::error:
            return "error";
        }
        throw std::invalid_argument{"unknown diagnostic severity"};
    };
    Json result = Json::array();
    for (const auto& diagnostic : diagnostics) {
        auto encoded = Json{{"severity", severity_name(diagnostic.severity)},
                            {"code", diagnostic.code},
                            {"message", diagnostic.message}};
        const auto add_index = [&encoded](const std::string_view name,
                                          const std::optional<std::size_t> index) {
            if (index.has_value()) {
                encoded[name] = *index;
            }
        };
        const auto& location = diagnostic.location;
        add_index("molecule_index", diagnostic.molecule_index);
        add_index("atom_index", location.atom_index);
        add_index("bond_index", location.bond_index);
        add_index("conformer_index", location.conformer_index);
        add_index("line", diagnostic.line);
        if (!location.atom_indices.empty()) {
            encoded["atom_indices"] = location.atom_indices;
        }
        if (!location.bond_indices.empty()) {
            encoded["bond_indices"] = location.bond_indices;
        }
        result.push_back(std::move(encoded));
    }
    return result;
}

[[nodiscard]] auto rounded(const double value, const double scale) -> double {
    return std::round(value * scale) / scale;
}

[[nodiscard]] auto import_metadata_json(const MoleculeImportMetadata& metadata) -> Json {
    auto policy = Json::object();
    const auto add_policy = [&policy](const std::string_view name,
                                      const std::optional<std::string>& value) {
        if (value.has_value()) {
            policy[name] = *value;
        }
    };
    add_policy("selection", metadata.policy.selection);
    add_policy("alternate_locations", metadata.policy.alternate_locations);
    add_policy("conformers", metadata.policy.conformers);
    add_policy("bonds", metadata.policy.bonds);

    return Json{{"format", to_string(metadata.format)},
                {"policy", std::move(policy)},
                {"source_connectivity", to_string(metadata.source_connectivity)}};
}

[[nodiscard]] auto portable_id_json(const PortableId& id) -> Json {
    if (id.empty()) {
        return nullptr;
    }
    return std::visit([](const auto& value) -> Json { return value; }, *id.value());
}

[[nodiscard]] auto portable_ids_json(const std::span<const PortableId> ids) -> Json {
    auto result = Json::array();
    for (const auto& id : ids) {
        result.push_back(portable_id_json(id));
    }
    return result;
}

[[nodiscard]] auto source_component_id(const ImportedMoleculeRecord& record,
                                       const std::size_t atom_index) -> std::optional<std::string> {
    if (!record.import_metadata.has_value()) {
        return std::nullopt;
    }

    const auto& metadata = *record.import_metadata;
    auto resolved = std::optional<std::string>{};
    auto ambiguous = false;
    const auto consider = [&resolved, &ambiguous](const SourceAtomReference& reference) {
        const auto value = detail::canonical_component_name(reference.structural_labels.get(), {});
        if (value.empty()) {
            return;
        }
        if (resolved.has_value() && *resolved != value) {
            ambiguous = true;
            return;
        }
        resolved = value;
    };

    consider(metadata.atoms.at(atom_index));
    for (const auto& conformer : metadata.conformers) {
        consider(conformer.sites.at(atom_index));
    }
    return ambiguous ? std::nullopt : resolved;
}

[[nodiscard]] auto fixed_ions_json(const calculation::FixedIons& fixed_ions,
                                   const std::span<const ImportedMoleculeRecord> records) -> Json {
    struct IonSummary {
        std::optional<std::string> id;
        double charge = 0.0;
        Json instances = Json::array();
    };
    auto groups = std::vector<IonSummary>{};
    for (const auto& source : fixed_ions.sources) {
        const auto id = source_component_id(records[source.molecule_index], source.atom_index);
        auto group = std::ranges::find_if(groups, [&id, &source](const auto& candidate) {
            return candidate.id == id && candidate.charge == source.charge;
        });
        if (group == groups.end()) {
            groups.push_back(IonSummary{.id = id, .charge = source.charge});
            group = std::prev(groups.end());
        }
        group->instances.push_back(
            {{"molecule_index", source.molecule_index}, {"atom_index", source.atom_index}});
    }

    auto result = Json::array();
    for (auto& group : groups) {
        result.push_back({{"component_id", group.id.has_value() ? Json(*group.id) : Json(nullptr)},
                          {"charge", group.charge},
                          {"instances", std::move(group.instances)}});
    }
    return result;
}

[[nodiscard]] auto record_json(const ImportedMoleculeRecord& record,
                               const ChargeCalculationResult& calculation_result,
                               const std::span<const charges::ChargeAssignment> record_charges,
                               const std::vector<ResultDiagnostic>& diagnostics) -> Json {
    Json input{{"source", record.identity.source}, {"record_index", record.identity.record_index}};
    if (!record.identity.record_id.empty()) {
        input["record_id"] = portable_id_json(record.identity.record_id);
    }
    if (record.caller_atom_ids.has_value()) {
        input["atom_ids"] = portable_ids_json(*record.caller_atom_ids);
    }
    if (record.import_metadata.has_value()) {
        input["import"] = import_metadata_json(*record.import_metadata);
    }

    Json result{{"input", std::move(input)}, {"status", to_string(calculation_result.status())}};
    if (calculation_result.charges() == nullptr) {
        result["diagnostics"] = diagnostics_json(diagnostics);
        return result;
    }

    Json assignments = Json::array();
    for (const auto& assignment : record_charges) {
        auto total_charge = 0.0;
        for (const auto value : assignment.charges.values()) {
            total_charge += value;
        }
        Json target{{"molecule_index", assignment.target.molecule_index}};
        if (assignment.target.conformer_index.has_value()) {
            target["conformer_index"] = *assignment.target.conformer_index;
        }
        Json encoded_assignment{
            {"scope", assignment.target.conformer_index.has_value() ? "conformer" : "molecule"},
            {"target", std::move(target)},
            {"charge_unit", "e"},
            {"charges", assignment.charges.values()},
            {"total_charge", total_charge}};
        assignments.push_back(std::move(encoded_assignment));
    }
    result["assignments"] = std::move(assignments);
    result["diagnostics"] = diagnostics_json(diagnostics);
    return result;
}

[[nodiscard]] auto provenance_json(const ChargeCalculationResult& result,
                                   const std::optional<ExecutionMetrics>& execution_metrics)
    -> Json {
    const auto optional_id = [](const std::optional<std::string>& id) -> Json {
        return id.has_value() ? Json{{"id", *id}} : Json(nullptr);
    };
    const auto option_value = [](const methods::MethodOptionValue& value) -> Json {
        return std::visit([](const auto& item) -> Json { return item; }, value);
    };
    const auto method_options_json =
        [&option_value](const std::map<std::string, methods::MethodOptions>& options) -> Json {
        Json encoded = Json::object();
        for (const auto& [method_id, values] : options) {
            Json method = Json::object();
            for (const auto& [id, value] : values.values()) {
                method[id] = option_value(value);
            }
            encoded[method_id] = std::move(method);
        }
        return encoded;
    };
    const auto threshold_value = [](const std::optional<std::size_t>& threshold) -> Json {
        return threshold.has_value() ? Json(*threshold) : Json("unlimited");
    };
    const auto optional_radius = [](const std::optional<double>& radius) -> Json {
        return radius.has_value() ? Json(*radius) : Json(nullptr);
    };

    const auto& requested_provenance = result.requested();
    Json requested{
        {"method", optional_id(requested_provenance.method_id)},
        {"parameter_set", optional_id(requested_provenance.parameter_set_id)},
        {"method_scope", calculation::to_string(requested_provenance.method_scope)},
        {"classification", {{"permissive_types", requested_provenance.permissive_types}}},
        {"resource_policy",
         {{"cutoff_atom_threshold", threshold_value(requested_provenance.cutoff_atom_threshold)},
          {"cover_atom_threshold", threshold_value(requested_provenance.cover_atom_threshold)},
          {"max_threads", requested_provenance.max_threads}}},
        {"execution",
         {{"kind", calculation::to_string(requested_provenance.execution_kind)},
          {"radius_angstrom", optional_radius(requested_provenance.execution_radius)}}}};
    requested["method_options"] = method_options_json(requested_provenance.method_options);

    auto effective = Json{{"method", nullptr},
                          {"parameter_set", nullptr},
                          {"warnings", Json::array()},
                          {"method_options", Json::object()}};
    if (const auto* execution = result.execution();
        execution != nullptr && execution->effective.has_value()) {
        const auto& value = *execution->effective;
        effective["method"] = optional_id(value.method_id);
        effective["parameter_set"] = optional_id(value.parameter_set_id);
        for (const auto& issue : value.warnings) {
            effective["warnings"].push_back(issue.message);
        }
        effective["execution"] = {
            {"mode", calculation::to_string(value.execution_policy.mode())},
            {"radius_angstrom", optional_radius(value.execution_policy.radius())}};
        effective["method_options"] =
            method_options_json({{value.method_id, value.method_options}});
        if (value.fixed_ions.has_value()) {
            effective["fixed_ions"] = fixed_ions_json(*value.fixed_ions, result.inputs());
        }
    }
    Json encoded{{"requested", std::move(requested)}, {"effective", std::move(effective)}};
    if (execution_metrics.has_value()) {
        const auto& metrics = *execution_metrics;
        encoded["execution_metrics"] = {
            {"started_at", metrics.started_at},
            {"ended_at", metrics.ended_at},
            {"runtime_seconds", rounded(metrics.runtime_seconds, metric_scale)},
            {"phases",
             {{"parsing_seconds", rounded(metrics.parsing_seconds, metric_scale)},
              {"applicability_seconds", rounded(metrics.applicability_seconds, metric_scale)},
              {"computation_seconds", rounded(metrics.computation_seconds, metric_scale)}}},
            {"peak_resident_memory_mb", rounded(metrics.peak_resident_memory_mb, metric_scale)}};
    }
    return encoded;
}

} // namespace

JsonWriter::JsonWriter(std::ostream& output) : output_{std::addressof(output)} {}

auto JsonWriter::write(const ChargeCalculationResult& result,
                       const std::optional<ExecutionMetrics>& execution_metrics) const -> void {
    Json records = Json::array();
    const auto inputs = result.inputs();
    for (std::size_t index = 0; index < inputs.size(); ++index) {
        records.push_back(record_json(inputs[index], result, record_assignments(result, index),
                                      charge_record_diagnostics(result, index)));
    }

    Json document{{"schema_version", "1.0"},
                  {"generator", {{"name", "ChargeFW"}, {"version", CHARGEFW_VERSION_STRING}}},
                  {"status", to_string(result.status())},
                  {"diagnostics", diagnostics_json(charge_result_diagnostics(result))},
                  {"results", std::move(records)}};
    document["calculation_provenance"] = provenance_json(result, execution_metrics);
    std::string serialized;
    try {
        serialized = document.dump(2);
    } catch (const Json::type_error& error) {
        // Error 316 is invalid UTF-8, such as a file name in another encoding on Linux.
        if (error.id != 316) {
            throw;
        }
        throw std::invalid_argument{"result JSON cannot hold text that is not valid UTF-8, such "
                                    "as a file name in another encoding"};
    }
    std::print(*output_, "{}\n", serialized);
    if (!*output_) {
        throw std::runtime_error{"failed to write JSON output"};
    }
}

} // namespace chargefw::adapters::native::json_output
