#pragma once

#include <chargefw/core/molecule.h>
#include <chargefw/methods/method_prerequisites.h>

#include "core/diagnostic_description.h"

#include <cstddef>
#include <span>
#include <string>
#include <string_view>

namespace chargefw::methods::detail {

[[nodiscard]] inline auto with_molecule_prefix(const core::Molecule& molecule,
                                               const std::size_t molecule_index,
                                               const PrerequisiteIssue& issue)
    -> PrerequisiteIssue {
    return PrerequisiteIssue{.kind = issue.kind,
                             .message =
                                 core::detail::molecule_description(molecule, molecule_index) +
                                 ": " + issue.message,
                             .molecule_index = molecule_index,
                             .location = issue.location};
}

// Why a parameter set does not apply to the atoms, which are in atom order and not empty. Fixed-ion
// planning rebuilds it for the original molecule, so it depends only on the atoms and the molecule.
[[nodiscard]] inline auto missing_atom_parameters_message(const core::Molecule& molecule,
                                                          const std::string_view parameter_set_id,
                                                          const std::span<const std::size_t> atoms)
    -> std::string {
    const auto parameter_set = "parameter set '" + std::string{parameter_set_id} + "'";
    if (atoms.size() == 1) {
        return parameter_set + " has no atom parameter matching " +
               core::detail::atom_description(molecule, atoms.front());
    }
    return parameter_set + " has no atom parameters for " + std::to_string(atoms.size()) + " of " +
           std::to_string(molecule.atom_count()) + " atoms (" +
           core::detail::element_counts_description(molecule, atoms) + "), such as " +
           core::detail::atom_description(molecule, atoms.front());
}

// The same for bonds, which are in bond order and not empty.
[[nodiscard]] inline auto missing_bond_parameters_message(const core::Molecule& molecule,
                                                          const std::string_view parameter_set_id,
                                                          const std::span<const std::size_t> bonds)
    -> std::string {
    const auto parameter_set = "parameter set '" + std::string{parameter_set_id} + "'";
    if (bonds.size() == 1) {
        return parameter_set + " has no bond parameter matching " +
               core::detail::bond_description(molecule, bonds.front());
    }
    return parameter_set + " has no bond parameters for " + std::to_string(bonds.size()) + " of " +
           std::to_string(molecule.bond_count()) + " bonds, such as " +
           core::detail::bond_description(molecule, bonds.front());
}

} // namespace chargefw::methods::detail
