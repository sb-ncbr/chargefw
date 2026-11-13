#include <chargefw/features/topology_features.h>

#include <stdexcept>

namespace chargefw::features {
namespace {

auto validate_atom_index(const core::Molecule& molecule, const std::size_t atom_index) -> void {
    if (atom_index >= molecule.atom_count()) {
        throw std::out_of_range{"atom index is outside the molecule"};
    }
}

} // namespace

TopologyFeatures::TopologyFeatures(const core::Molecule& molecule)
    : molecule_{&molecule}, neighbor_indices_(molecule.atom_count()),
      incident_bond_indices_(molecule.atom_count()) {
    const auto bonds = molecule.bonds();

    for (std::size_t bond_index = 0; bond_index < bonds.size(); ++bond_index) {
        const auto& bond = bonds[bond_index];

        const auto first_atom_index = bond.first_atom_index();
        const auto second_atom_index = bond.second_atom_index();

        neighbor_indices_[first_atom_index].push_back(second_atom_index);
        incident_bond_indices_[first_atom_index].push_back(bond_index);

        neighbor_indices_[second_atom_index].push_back(first_atom_index);
        incident_bond_indices_[second_atom_index].push_back(bond_index);
    }
}

auto TopologyFeatures::molecule() const noexcept -> const core::Molecule& {
    return *molecule_;
}

auto TopologyFeatures::degree(const std::size_t atom_index) const -> std::size_t {
    validate_atom_index(molecule(), atom_index);
    return neighbor_indices_[atom_index].size();
}

auto TopologyFeatures::neighbor_indices(const std::size_t atom_index) const
    -> std::span<const std::size_t> {
    validate_atom_index(molecule(), atom_index);

    const auto& neighbors = neighbor_indices_[atom_index];
    return {neighbors.data(), neighbors.size()};
}

auto TopologyFeatures::incident_bond_indices(const std::size_t atom_index) const
    -> std::span<const std::size_t> {
    validate_atom_index(molecule(), atom_index);

    const auto& bonds = incident_bond_indices_[atom_index];
    return {bonds.data(), bonds.size()};
}

} // namespace chargefw::features
