#include <chargefw/core/periodic_table.h>
#include <chargefw/parameters/classification/parameter_classifier.h>

#include "core/diagnostic_description.h"
#include "parameters/classification/classification_checks.h"
#include "parameters/classification/type_keys.h"

#include <algorithm>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace chargefw::parameters {
namespace {

auto permissive_bond_order_type(const core::BondOrder order) -> std::string {
    const auto value = static_cast<int>(core::bond_order_value(order));
    return std::to_string(std::max(1, value - 1));
}

auto highest_bond_order_rank(const features::TopologyFeatures& topology,
                             const std::size_t atom_index) -> int {
    const auto& molecule = topology.molecule();
    auto highest = 0;

    for (const auto bond_index : topology.incident_bond_indices(atom_index)) {
        highest = std::max(
            highest, static_cast<int>(core::bond_order_value(molecule.bond(bond_index).order())));
    }

    return highest;
}

auto highest_bond_order_type(const features::TopologyFeatures& topology,
                             const std::size_t atom_index, const bool permissive) -> std::string {
    const auto highest = highest_bond_order_rank(topology, atom_index);

    // An atom without bonds, such as the oxygen of water without its hydrogens, falls back to
    // single-bonded parameters.
    if (permissive) {
        return std::to_string(std::max(1, highest - 1));
    }

    return std::to_string(highest);
}

auto bonded_elements_type(const features::TopologyFeatures& topology, const std::size_t atom_index)
    -> std::string {
    const auto& molecule = topology.molecule();
    std::vector<std::string_view> symbols;
    symbols.reserve(topology.degree(atom_index));

    for (const auto neighbor_index : topology.neighbor_indices(atom_index)) {
        symbols.push_back(core::element_symbol(molecule.atom(neighbor_index).atomic_number()));
    }

    std::ranges::sort(symbols);

    std::string type;

    for (const auto symbol : symbols) {
        type += symbol;
    }

    return type;
}

auto atom_type_for(const features::TopologyFeatures& topology, const std::size_t atom_index,
                   const AtomParameterClassificationKind classification, const bool permissive)
    -> std::string {
    switch (classification) {
    case AtomParameterClassificationKind::plain:
        return "*";

    case AtomParameterClassificationKind::highest_bond_order:
        return highest_bond_order_type(topology, atom_index, permissive);

    case AtomParameterClassificationKind::bonded_elements:
        return bonded_elements_type(topology, atom_index);
    }

    throw std::logic_error{"unknown atom parameter classification"};
}

auto bond_type_for(const core::Molecule& molecule, const std::size_t bond_index,
                   const BondParameterClassificationKind classification, const bool permissive)
    -> std::string {
    switch (classification) {
    case BondParameterClassificationKind::plain:
        return "*";

    case BondParameterClassificationKind::bond_order:
        if (permissive) {
            return permissive_bond_order_type(molecule.bond(bond_index).order());
        }

        return std::to_string(core::bond_order_value(molecule.bond(bond_index).order()));
    }

    throw std::logic_error{"unknown bond parameter classification"};
}

auto matches_atom_key(const features::TopologyFeatures& topology, const std::size_t atom_index,
                      const AtomParameterKey& key, const bool permissive) -> bool {
    const auto& atom = topology.molecule().atom(atom_index);

    if (key.atomic_number != 0 && key.atomic_number != atom.atomic_number()) {
        return false;
    }

    return key.type == atom_type_for(topology, atom_index, key.classification, permissive);
}

auto matches_bond_key(const features::TopologyFeatures& topology, const std::size_t bond_index,
                      const BondParameterKey& key, const bool permissive) -> bool {
    const auto& molecule = topology.molecule();
    const auto& bond = molecule.bond(bond_index);

    const auto first_atom_index = bond.first_atom_index();
    const auto second_atom_index = bond.second_atom_index();

    const auto forward = matches_atom_key(topology, first_atom_index, key.first_atom, permissive) &&
                         matches_atom_key(topology, second_atom_index, key.second_atom, permissive);

    const auto reverse =
        matches_atom_key(topology, first_atom_index, key.second_atom, permissive) &&
        matches_atom_key(topology, second_atom_index, key.first_atom, permissive);

    if (!forward && !reverse) {
        return false;
    }

    return key.bond.type ==
           bond_type_for(molecule, bond_index, key.bond.classification, permissive);
}

template <typename Parameters, typename Matches>
auto find_entry(const Parameters& parameters, const ClassificationOptions& options, Matches matches)
    -> std::optional<std::size_t> {
    const auto scan = [&](const bool permissive) -> std::optional<std::size_t> {
        for (std::size_t entry_index = 0; entry_index < parameters.size(); ++entry_index) {
            if (matches(parameters[entry_index].key, permissive)) {
                return entry_index;
            }
        }
        return std::nullopt;
    };
    auto match = scan(false);

    if (!match && options.permissive_types) {
        match = scan(true);
    }

    return match;
}

} // namespace

namespace detail {

auto atom_key_less(const AtomParameterKey& first, const AtomParameterKey& second) -> bool {
    return std::tie(first.atomic_number, first.classification, first.type) <
           std::tie(second.atomic_number, second.classification, second.type);
}

auto atom_type_key(const features::TopologyFeatures& topology, const std::size_t atom_index,
                   const AtomParameterClassificationKind classification) -> AtomParameterKey {
    return {.atomic_number = topology.molecule().atom(atom_index).atomic_number(),
            .classification = classification,
            .type = atom_type_for(topology, atom_index, classification, false)};
}

auto bond_type_key(const features::TopologyFeatures& topology, const std::size_t bond_index,
                   const AtomParameterClassificationKind atom_classification,
                   const BondParameterClassificationKind bond_classification) -> BondParameterKey {
    const auto& molecule = topology.molecule();
    const auto& bond = molecule.bond(bond_index);
    auto lower = atom_type_key(topology, bond.first_atom_index(), atom_classification);
    auto upper = atom_type_key(topology, bond.second_atom_index(), atom_classification);
    if (atom_key_less(upper, lower)) {
        std::swap(lower, upper);
    }
    return {.first_atom = std::move(lower),
            .second_atom = std::move(upper),
            .bond = {.classification = bond_classification,
                     .type = bond_type_for(molecule, bond_index, bond_classification, false)}};
}

} // namespace detail

auto try_classify_parameters(const features::TopologyFeatures& topology,
                             const ParameterSet& parameters, const ClassificationOptions& options)
    -> ClassificationResult {
    const auto& molecule = topology.molecule();
    std::vector<ClassificationIssue> issues;
    std::vector<std::size_t> atom_indices;
    std::vector<std::size_t> bond_indices;

    if (!parameters.atom().empty()) {
        atom_indices.reserve(molecule.atom_count());

        for (std::size_t atom_index = 0; atom_index < molecule.atom_count(); ++atom_index) {
            const auto match =
                find_entry(parameters.atom(), options,
                           [&](const AtomParameterKey& key, const bool permissive) {
                               return matches_atom_key(topology, atom_index, key, permissive);
                           });

            if (!match) {
                issues.push_back(ClassificationIssue{
                    .kind = ClassificationIssueKind::missing_atom_parameter,
                    .object_index = atom_index,
                    .message = "parameter set '" + std::string{parameters.id()} +
                               "' has no atom parameter matching " +
                               core::detail::atom_description(molecule, atom_index)});
                continue;
            }

            atom_indices.push_back(*match);
        }
    }

    if (!parameters.bond().empty()) {
        bond_indices.reserve(molecule.bond_count());

        for (std::size_t bond_index = 0; bond_index < molecule.bond_count(); ++bond_index) {
            const auto match =
                find_entry(parameters.bond(), options,
                           [&](const BondParameterKey& key, const bool permissive) {
                               return matches_bond_key(topology, bond_index, key, permissive);
                           });

            if (!match) {
                issues.push_back(ClassificationIssue{
                    .kind = ClassificationIssueKind::missing_bond_parameter,
                    .object_index = bond_index,
                    .message = "parameter set '" + std::string{parameters.id()} +
                               "' has no bond parameter matching " +
                               core::detail::bond_description(molecule, bond_index)});
                continue;
            }

            bond_indices.push_back(*match);
        }
    }

    if (!issues.empty()) {
        return ClassificationResult{std::move(issues)};
    }

    auto classification =
        ParameterClassification{AtomParameterClassification{std::move(atom_indices)},
                                BondParameterClassification{std::move(bond_indices)}};

    validate_parameter_classification(molecule, parameters, classification);

    return ClassificationResult{std::move(classification)};
}

auto classify_parameters(const features::TopologyFeatures& topology, const ParameterSet& parameters,
                         const ClassificationOptions& options) -> ParameterClassification {
    const auto result = try_classify_parameters(topology, parameters, options);

    if (!result) {
        throw std::invalid_argument{result.issues().front().message};
    }

    return result.classification();
}

} // namespace chargefw::parameters
