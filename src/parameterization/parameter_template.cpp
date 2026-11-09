#include <chargefw/parameterization/parameter_template.h>

#include "parameters/classification/type_keys.h"

#include <chargefw/parameters/models/atom_parameters.h>
#include <chargefw/parameters/models/bond_parameters.h>
#include <chargefw/parameters/models/common_parameters.h>

#include <cstddef>
#include <set>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace chargefw::parameterization {
namespace {

using parameters::AtomParameterKey;
using parameters::BondParameterKey;

struct AtomKeyLess {
    [[nodiscard]] auto operator()(const AtomParameterKey& first,
                                  const AtomParameterKey& second) const -> bool {
        return parameters::detail::atom_key_less(first, second);
    }
};

// Orders bond keys by their first atom key, then their second, then classifier and type.
struct BondKeyLess {
    [[nodiscard]] auto operator()(const BondParameterKey& first,
                                  const BondParameterKey& second) const -> bool {
        const auto atom_less = AtomKeyLess{};
        if (first.first_atom != second.first_atom) {
            return atom_less(first.first_atom, second.first_atom);
        }
        if (first.second_atom != second.second_atom) {
            return atom_less(first.second_atom, second.second_atom);
        }
        return std::tie(first.bond.classification, first.bond.type) <
               std::tie(second.bond.classification, second.bond.type);
    }
};

auto zero_parameters(const std::vector<std::string>& names)
    -> std::vector<parameters::NamedParameter> {
    auto result = std::vector<parameters::NamedParameter>{};
    result.reserve(names.size());
    for (const auto& name : names) {
        result.push_back({.name = name, .value = 0.0});
    }
    return result;
}

} // namespace

auto make_parameter_template(const features::PreparedMoleculeCollection& molecules,
                             const ParameterNames& names,
                             const parameters::AtomParameterClassificationKind atom_types,
                             const parameters::BondParameterClassificationKind bond_types,
                             parameters::ParameterSetMetadata metadata)
    -> parameters::ParameterSet {
    auto atom_keys = std::set<AtomParameterKey, AtomKeyLess>{};
    auto bond_keys = std::set<BondParameterKey, BondKeyLess>{};
    for (const auto& prepared : molecules.molecules()) {
        const auto& molecule = prepared.molecule();
        const auto& topology = prepared.topology();
        if (!names.atom.empty()) {
            for (std::size_t atom_index = 0; atom_index < molecule.atom_count(); ++atom_index) {
                atom_keys.insert(
                    parameters::detail::atom_type_key(molecule, topology, atom_index, atom_types));
            }
        }
        if (!names.bond.empty()) {
            for (std::size_t bond_index = 0; bond_index < molecule.bond_count(); ++bond_index) {
                bond_keys.insert(parameters::detail::bond_type_key(molecule, topology, bond_index,
                                                                   atom_types, bond_types));
            }
        }
    }

    auto atom_entries = std::vector<parameters::AtomParameterEntry>{};
    atom_entries.reserve(atom_keys.size());
    for (const auto& key : atom_keys) {
        atom_entries.push_back({.key = key, .parameters = zero_parameters(names.atom)});
    }
    auto bond_entries = std::vector<parameters::BondParameterEntry>{};
    bond_entries.reserve(bond_keys.size());
    for (const auto& key : bond_keys) {
        bond_entries.push_back({.key = key, .parameters = zero_parameters(names.bond)});
    }
    return parameters::ParameterSet{std::move(metadata),
                                    parameters::CommonParameters{zero_parameters(names.common)},
                                    parameters::AtomParameters{std::move(atom_entries)},
                                    parameters::BondParameters{std::move(bond_entries)}};
}

} // namespace chargefw::parameterization
