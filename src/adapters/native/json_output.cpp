#include <chargefw/adapters/native/json_output.h>

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
        add_index("molecule_index", diagnostic.molecule_index);
        add_index("atom_index", diagnostic.atom_index);
        add_index("bond_index", diagnostic.bond_index);
        add_index("conformer_index", diagnostic.conformer_index);
        add_index("line", diagnostic.line);
        result.push_back(std::move(encoded));
    }
    return result;
}

[[nodiscard]] auto rounded(const double value, const double scale) -> double {
    return std::round(value * scale) / scale;
}

[[nodiscard]] auto source_format_name(const MolecularSourceFormat format) -> std::string_view {
    switch (format) {
    case MolecularSourceFormat::molecule_json:
        return "molecule-json";
    case MolecularSourceFormat::mol:
        return "mol";
    case MolecularSourceFormat::sdf:
        return "sdf";
    case MolecularSourceFormat::mol2:
        return "mol2";
    case MolecularSourceFormat::pdb:
        return "pdb";
    case MolecularSourceFormat::mmcif:
        return "mmcif";
    }
    throw std::invalid_argument{"unknown molecular source format"};
}

[[nodiscard]] auto connectivity_name(const SourceConnectivity connectivity) -> std::string_view {
    switch (connectivity) {
    case SourceConnectivity::absent:
        return "absent";
    case SourceConnectivity::explicitly_empty:
        return "explicitly-empty";
    case SourceConnectivity::present:
        return "present";
    }
    throw std::invalid_argument{"unknown source connectivity state"};
}

[[nodiscard]] auto import_metadata_json(const MoleculeImportMetadata& metadata) -> Json {
    auto policy = Json::object();
    const auto add_policy = [&policy](const std::string_view name,
                                      const std::optional<std::string>& value) {
        if (value.has_value()) {
            policy[name] = *value;
        }
    };
    add_policy("record_selection", metadata.record_selection);
    add_policy("alternate_location_selection", metadata.alternate_location_selection);
    add_policy("conformer_selection", metadata.conformer_selection);
    add_policy("bond_strategy", metadata.bond_strategy);

    return Json{{"format", source_format_name(metadata.format)},
                {"policy", std::move(policy)},
                {"source_connectivity", connectivity_name(metadata.source_connectivity)}};
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

[[nodiscard]] auto component_id(const SourceAtomReference& reference)
    -> std::optional<std::string> {
    if (!reference.structural_labels.has_value()) {
        return std::nullopt;
    }
    const auto& labels = *reference.structural_labels;
    if (labels.label.residue.has_value() && !labels.label.residue->empty()) {
        return labels.label.residue;
    }
    if (labels.author.residue.has_value() && !labels.author.residue->empty()) {
        return labels.author.residue;
    }
    return std::nullopt;
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
        const auto value = component_id(reference);
        if (!value.has_value()) {
            return;
        }
        if (resolved.has_value() && *resolved != *value) {
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
                               const calculation::ExecutionResult& execution,
                               const std::size_t molecule_index,
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

    Json result{{"input", std::move(input)}, {"status", calculation::to_string(execution.status)}};
    if (!execution.charges.has_value()) {
        result["diagnostics"] = diagnostics_json(diagnostics);
        return result;
    }

    Json assignments = Json::array();
    for (const auto& assignment : execution.charges->assignments()) {
        if (assignment.target.molecule_index != molecule_index) {
            continue;
        }
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
        {"classification", {{"permissive_types", requested_provenance.permissive_types}}},
        {"resource_policy",
         {{"cutoff_atom_threshold", threshold_value(requested_provenance.cutoff_atom_threshold)},
          {"cover_atom_threshold", threshold_value(requested_provenance.cover_atom_threshold)},
          {"max_threads", requested_provenance.max_threads}}},
        {"execution",
         {{"kind", requested_provenance.execution_kind},
          {"radius_angstrom", optional_radius(requested_provenance.execution_radius)}}}};
    requested["method_options"] = method_options_json(requested_provenance.method_options);

    auto effective = Json{{"method", nullptr},
                          {"parameter_set", nullptr},
                          {"warnings", Json::array()},
                          {"method_options", Json::object()}};
    if (result.execution().effective.has_value()) {
        const auto& value = *result.execution().effective;
        effective["method"] = optional_id(value.method_id);
        effective["parameter_set"] = optional_id(value.parameter_set_id);
        for (const auto& issue : value.execution_issues) {
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
                       const std::string_view generator_version,
                       const std::optional<ExecutionMetrics>& execution_metrics) const -> void {
    Json records = Json::array();
    const auto inputs = result.inputs();
    for (std::size_t index = 0; index < inputs.size(); ++index) {
        records.push_back(record_json(inputs[index], result.execution(), index,
                                      charge_record_diagnostics(result, index)));
    }

    Json document{{"schema_version", "1.0"},
                  {"generator", {{"name", "ChargeFW"}, {"version", generator_version}}},
                  {"status", calculation::to_string(result.execution().status)},
                  {"diagnostics", diagnostics_json(charge_result_diagnostics(result))},
                  {"results", std::move(records)}};
    document["calculation_provenance"] = provenance_json(result, execution_metrics);
    std::print(*output_, "{}\n", document.dump(2));
}

} // namespace chargefw::adapters::native::json_output
