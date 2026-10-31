#include "cli_support.h"

#include <chargefw/core/periodic_table.h>
#include <chargefw/methods/method_registry.h>
#include <chargefw/parameters/io/parameter_set_io.h>
#include <chargefw/parameters/models/parameter_key.h>

#include <algorithm>
#include <array>
#include <format>
#include <iostream>
#include <map>
#include <print>
#include <stdexcept>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

namespace chargefw::cli {
namespace {

struct SummaryRow {
    std::string_view id;
    std::string_view name;
};

void print_summary(const std::vector<SummaryRow>& rows) {
    std::size_t id_width = 0;
    for (const auto& row : rows) {
        id_width = std::max(id_width, row.id.size());
    }
    for (const auto& row : rows) {
        std::println("{:<{}}  {}", row.id, id_width, row.name);
    }
}

[[nodiscard]] auto format_bytes(double bytes) -> std::string {
    constexpr auto units = std::array{"B", "KiB", "MiB", "GiB", "TiB"};
    auto unit = std::size_t{0};
    while (bytes >= 1024.0 && unit + 1 < units.size()) {
        bytes /= 1024.0;
        ++unit;
    }
    return std::format("{:.1f}{}", bytes, units[unit]);
}

// Shows a classifier and its expected type; an empty type is quoted to stay visible.
[[nodiscard]] auto classified_type(const std::string_view classifier, const std::string& type)
    -> std::string {
    return std::format("{}={}", classifier, type.empty() ? std::string{"\"\""} : type);
}

struct AtomKeyText {
    std::string element;
    std::string type;
};

[[nodiscard]] auto atom_key_text(const parameters::AtomParameterKey& key) -> AtomKeyText {
    return {.element = key.atomic_number == 0
                           ? std::string{"*"}
                           : std::string{core::element_symbol(key.atomic_number)},
            .type = classified_type(parameters::to_string(key.classification), key.type)};
}

// Classifier names in order of first use, joined by commas.
template <typename Entries, typename Classifier>
[[nodiscard]] auto classifier_names(const Entries& entries, const Classifier& classifier)
    -> std::string {
    auto names = std::vector<std::string_view>{};
    for (const auto& entry : entries) {
        const auto name = classifier(entry);
        if (std::ranges::find(names, name) == names.end()) {
            names.push_back(name);
        }
    }
    auto joined = std::string{};
    for (const auto name : names) {
        joined += joined.empty() ? "" : ",";
        joined += name;
    }
    return joined;
}

void print_atom_keys(const parameters::AtomParameters& atoms, const bool keys) {
    std::print("atom keys: {}", atoms.size());
    if (!atoms.empty()) {
        std::print(" classifiers={}", classifier_names(atoms.entries(), [](const auto& entry) {
                       return parameters::to_string(entry.key.classification);
                   }));
    }
    std::println();
    if (!keys) {
        return;
    }
    auto rows = std::vector<AtomKeyText>{};
    auto element_width = std::size_t{0};
    for (const auto& entry : atoms.entries()) {
        rows.push_back(atom_key_text(entry.key));
        element_width = std::max(element_width, rows.back().element.size());
    }
    const auto index_width = std::to_string(rows.size()).size();
    for (std::size_t index = 0; index < rows.size(); ++index) {
        std::println("  {:>{}}  {:<{}}  {}", index + 1, index_width, rows[index].element,
                     element_width, rows[index].type);
    }
}

void print_bond_keys(const parameters::BondParameters& bonds, const bool keys) {
    std::print("bond keys: {}", bonds.size());
    if (!bonds.empty()) {
        std::print(" atom-classifiers={} bond-classifiers={}",
                   classifier_names(bonds.entries(),
                                    [](const auto& entry) {
                                        return parameters::to_string(
                                            entry.key.first_atom.classification);
                                    }),
                   classifier_names(bonds.entries(), [](const auto& entry) {
                       return parameters::to_string(entry.key.bond.classification);
                   }));
    }
    std::println();
    if (!keys) {
        return;
    }
    struct Row {
        AtomKeyText first;
        AtomKeyText second;
        std::string bond;
    };
    auto rows = std::vector<Row>{};
    auto element_width = std::size_t{0};
    auto type_width = std::size_t{0};
    for (const auto& entry : bonds.entries()) {
        rows.push_back(
            {.first = atom_key_text(entry.key.first_atom),
             .second = atom_key_text(entry.key.second_atom),
             .bond = classified_type(parameters::to_string(entry.key.bond.classification),
                                     entry.key.bond.type)});
        for (const auto* atom : {&rows.back().first, &rows.back().second}) {
            element_width = std::max(element_width, atom->element.size());
            type_width = std::max(type_width, atom->type.size());
        }
    }
    const auto index_width = std::to_string(rows.size()).size();
    for (std::size_t index = 0; index < rows.size(); ++index) {
        const auto& row = rows[index];
        std::println("  {:>{}}  {:<{}}  {:<{}}  {:<{}}  {:<{}}  {}", index + 1, index_width,
                     row.first.element, element_width, row.first.type, type_width,
                     row.second.element, element_width, row.second.type, type_width, row.bond);
    }
}

} // namespace

void print_inspection(const ImportedCollection& imported) {
    std::println("records: {}", imported.molecules.size());
    for (std::size_t index = 0; index < imported.molecules.size(); ++index) {
        const auto& molecule = imported.molecules[index];
        std::map<int, std::size_t> elements;
        for (const auto& atom : molecule.atoms()) {
            ++elements[atom.atomic_number()];
        }
        std::print("record {} ({}) atoms={} bonds={} conformers={} coordinates={} formal_charge={} "
                   "elements=",
                   index, imported.records[index].identity.record_id.display_string(),
                   molecule.atom_count(), molecule.bond_count(), molecule.conformer_count(),
                   molecule.has_coordinates(), core::total_formal_charge(molecule));
        bool first = true;
        for (const auto& [atomic_number, count] : elements) {
            std::print("{}{}:{}", first ? "" : ",", core::element_symbol(atomic_number), count);
            first = false;
        }
        std::println();
    }
}

void print_applicability(const calculation::AssessmentResult& assessment,
                         const std::size_t max_threads) {
    std::println("runnable plans: {}", assessment.plans().size());
    for (const auto& plan : assessment.plans()) {
        const auto estimate = calculation::estimate(assessment, plan, max_threads);
        std::println("plan method={} parameter_set={} execution={} memory={} operations={:.1e}",
                     plan.method().id(),
                     plan.parameter_set() == nullptr ? "-" : plan.parameter_set()->id(),
                     calculation::to_string(plan.policy().mode()),
                     format_bytes(estimate.memory_bytes), estimate.operations);
    }
    std::println("rejected alternatives: {}", assessment.rejections().size());
    for (const auto& rejected : assessment.rejections()) {
        std::print("rejected method={}", rejected.method_id);
        if (rejected.parameter_set_id.has_value()) {
            std::print(" parameter_set={}", *rejected.parameter_set_id);
        }
        if (rejected.policy.has_value()) {
            std::print(" execution={}", calculation::to_string(rejected.policy->mode()));
        }
        for (const auto& issue : rejected.issues) {
            std::visit([](const auto& value) { std::print("; {}", value.message); }, issue);
        }
        std::println();
    }
    const auto* plan = assessment.default_plan();
    if (plan == nullptr) {
        std::println("selected execution: none");
        return;
    }
    std::println("selected method={} parameter_set={} execution={}", plan->method().id(),
                 plan->parameter_set() == nullptr ? "-" : plan->parameter_set()->id(),
                 calculation::to_string(plan->policy().mode()));
}

void print_methods(const std::string& method_id) {
    const auto& registry = methods::method_registry();
    if (method_id.empty()) {
        auto rows = std::vector<SummaryRow>{};
        rows.reserve(registry.methods().size());
        for (const auto& method : registry.methods()) {
            rows.push_back({method->id(), method->metadata().full_name});
        }
        print_summary(rows);
        return;
    }

    const auto* method = registry.find(method_id);
    if (method == nullptr) {
        throw std::invalid_argument{"method '" + method_id + "' is not registered"};
    }

    const auto& metadata = method->metadata();
    const auto requirements = method->requirements();
    std::println("id: {}", metadata.id);
    std::println("name: {}", metadata.name);
    std::println("full name: {}", metadata.full_name);
    std::println("publication: {}", metadata.publication.value_or("-"));
    std::println("notes: {}", metadata.notes.empty() ? "-" : metadata.notes);
    std::println("recommended: {}", metadata.recommended ? "yes" : "no");
    std::println("requires coordinates: {}", requirements.coordinates ? "yes" : "no");
    std::println("requires bonds: {}", requirements.bonds ? "yes" : "no");
    std::println("time complexity: {}", methods::complexity_notation(requirements.resources.time));
    std::println("memory complexity: {}",
                 methods::complexity_notation(requirements.resources.memory));
    std::println("supports cutoff: {}", requirements.resources.supports_cutoff ? "yes" : "no");
    std::println("supports cover: {}", requirements.resources.supports_cover ? "yes" : "no");
    std::println("supports fixed ions: {}",
                 requirements.supports_fixed_point_sources ? "yes" : "no");
    std::println("options:{}", method->option_schema().empty() ? " none" : "");
    for (const auto& option : method->option_schema()) {
        std::print("  {} (default=", option.id);
        std::visit([](const auto& value) { std::print("{}", value); }, option.default_value);
        std::print(")");
        if (!option.choices.empty()) {
            std::print(" choices=");
            for (std::size_t index = 0; index < option.choices.size(); ++index) {
                if (index != 0) {
                    std::print(",");
                }
                std::visit([](const auto& value) { std::print("{}", value); },
                           option.choices[index]);
            }
        }
        if (option.minimum.has_value()) {
            std::print("{}", option.minimum_inclusive ? " minimum>=" : " minimum>");
            std::visit([](const auto& value) { std::print("{}", value); }, *option.minimum);
        }
        if (option.maximum.has_value()) {
            std::print("{}", option.maximum_inclusive ? " maximum<=" : " maximum<");
            std::visit([](const auto& value) { std::print("{}", value); }, *option.maximum);
        }
        std::println();
    }
}

void print_parameter_sets(const std::string& parameter_set_id, const std::string& method_id,
                          const bool keys) {
    const auto parameter_sets = parameters::load_default_parameter_sets();
    if (parameter_set_id.empty()) {
        if (!method_id.empty() && methods::method_registry().find(method_id) == nullptr) {
            throw std::invalid_argument{"method '" + method_id + "' is not registered"};
        }
        const auto visible = [&method_id](const parameters::ParameterSet& parameter_set) {
            return method_id.empty() || parameter_set.method_id() == method_id;
        };
        auto rows = std::vector<SummaryRow>{};
        rows.reserve(parameter_sets.size());
        for (const auto& parameter_set : parameter_sets) {
            if (visible(parameter_set)) {
                rows.push_back({parameter_set.id(), parameter_set.name()});
            }
        }
        print_summary(rows);
        return;
    }
    const auto parameter_set =
        std::ranges::find(parameter_sets, parameter_set_id, &parameters::ParameterSet::id);
    if (parameter_set == parameter_sets.end()) {
        throw std::invalid_argument{"parameter set '" + parameter_set_id + "' was not found"};
    }

    const auto& metadata = parameter_set->metadata();
    std::println("id: {}", metadata.id);
    std::println("method: {}", metadata.method_id);
    std::println("name: {}", metadata.name);
    std::println("publication: {}", metadata.publication.empty() ? "-" : metadata.publication);
    std::println("notes: {}", metadata.notes.empty() ? "-" : metadata.notes);
    print_atom_keys(parameter_set->atom(), keys);
    print_bond_keys(parameter_set->bond(), keys);
}

} // namespace chargefw::cli
