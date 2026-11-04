#include <chargefw/calculation/input_issues.h>

#include "core/diagnostic_description.h"
#include "input_checks.h"

#include <chargefw/features/prepared_molecule.h>

#include <algorithm>
#include <format>
#include <stdexcept>
#include <string>

namespace chargefw::calculation {

auto to_string(const InputIssueKind value) -> std::string_view {
    switch (value) {
    case InputIssueKind::hydrogens_missing:
        return "hydrogens_missing";
    case InputIssueKind::atoms_too_close:
        return "atoms_too_close";
    }
    throw std::invalid_argument{"unknown input issue kind"};
}

namespace detail {

auto find_input_issues(const core::MoleculeCollection& molecules) -> std::vector<InputIssue> {
    constexpr auto hydrogen = 1;
    constexpr auto carbon = 6;
    auto result = std::vector<InputIssue>{};
    for (std::size_t index = 0; index < molecules.size(); ++index) {
        const auto atoms = molecules[index].atoms();
        const auto has = [atoms](const int atomic_number) {
            return std::ranges::any_of(atoms, [atomic_number](const auto& atom) {
                return atom.atomic_number() == atomic_number;
            });
        };
        // Hydrogens are missing from a crystal structure before protonation or from an RDKit
        // molecule without AddHs(); a molecule such as CO2 or CCl4 has none, so this only warns.
        if (has(carbon) && !has(hydrogen)) {
            result.push_back(InputIssue{
                .kind = InputIssueKind::hydrogens_missing,
                .message = "the molecule contains carbon but no hydrogen atoms; if its hydrogens "
                           "were not added, add them before calculating",
                .molecule_index = index});
        }
    }
    return result;
}

auto add_close_atom_issues(const std::span<const InputIssue> assessment_warnings,
                           const core::MoleculeCollection& molecules,
                           const features::PreparedMoleculeCollection& prepared,
                           const FixedChargePartition* fixed_charge_partition)
    -> std::vector<InputIssue> {
    auto result = std::vector<InputIssue>{};
    auto next_warning = assessment_warnings.begin();
    for (std::size_t index = 0; index < molecules.size(); ++index) {
        for (; next_warning != assessment_warnings.end() && next_warning->molecule_index == index;
             ++next_warning) {
            result.push_back(*next_warning);
        }
        const auto& molecule = molecules[index];
        // Without fixed ions, the active part is a copy with the same atom indices.
        const auto searched_anew = fixed_charge_partition != nullptr &&
                                   !fixed_charge_partition->targets[index].sources.empty();
        for (std::size_t conformer = 0; conformer < molecule.conformer_count(); ++conformer) {
            const auto close =
                searched_anew
                    ? features::find_close_atoms(features::ConformerFeatures{molecule, conformer})
                    : prepared[index].geometry_defects()[conformer].close_atoms;
            if (!close.has_value()) {
                continue;
            }
            const auto [first, second] = close->closest;
            const auto pair =
                std::format("{} and {}", core::detail::atom_description(molecule, first),
                            core::detail::atom_description(molecule, second));
            const auto finding =
                close->pair_count == 1
                    ? std::format("{} are {:.2f} Å apart, closer than {} times the sum of their "
                                  "covalent radii",
                                  pair, close->distance, features::close_atom_factor)
                    : std::format("{} pairs of atoms are closer than {} times the sum of their "
                                  "covalent radii; closest relative to that sum are {}, {:.2f} Å "
                                  "apart",
                                  close->pair_count, features::close_atom_factor, pair,
                                  close->distance);
            result.push_back(InputIssue{
                .kind = InputIssueKind::atoms_too_close,
                .message = std::format("conformer {}: {}; check the coordinates, because charges "
                                       "calculated from them may be wrong",
                                       conformer + 1, finding),
                .molecule_index = index,
                .conformer_index = conformer,
                .atom_index = second,
                .atom_indices = close->atoms});
        }
    }
    return result;
}

} // namespace detail
} // namespace chargefw::calculation
