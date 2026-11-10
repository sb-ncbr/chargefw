#include <chargefw/parameterization/parameter_entries.h>

#include "calculation/assessment_access.h"

#include <chargefw/core/molecule.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/methods/method.h>
#include <chargefw/parameters/classification/parameter_classification.h>
#include <chargefw/parameters/models/parameter_set.h>

#include <algorithm>
#include <cstddef>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace chargefw::parameterization {
namespace {

using calculation::detail::AssessmentAccess;

[[nodiscard]] auto parameter_set(const calculation::ExecutionPlan& plan)
    -> const parameters::ParameterSet& {
    const auto* set = plan.parameter_set();
    if (set == nullptr) {
        throw std::invalid_argument{"method '" + std::string{plan.method().id()} +
                                    "' has no parameters"};
    }
    return *set;
}

} // namespace

auto entry_counts(const calculation::AssessmentResult& assessment,
                  const calculation::ExecutionPlan& plan) -> EntryCounts {
    const auto classifications = AssessmentAccess::classifications(assessment, plan);
    const auto& set = parameter_set(plan);
    auto result = EntryCounts{.atoms = std::vector<std::size_t>(set.atom().size(), 0),
                              .bonds = std::vector<std::size_t>(set.bond().size(), 0)};
    for (const auto& classification : classifications) {
        for (const auto entry : classification.atom().parameter_entry_indices()) {
            ++result.atoms[entry];
        }
        for (const auto entry : classification.bond().parameter_entry_indices()) {
            ++result.bonds[entry];
        }
    }
    return result;
}

auto atom_entries(const calculation::AssessmentResult& assessment,
                  const calculation::ExecutionPlan& plan)
    -> std::vector<std::vector<std::optional<std::size_t>>> {
    const auto classifications = AssessmentAccess::classifications(assessment, plan);
    static_cast<void>(parameter_set(plan));
    const auto& molecules = assessment.molecules();
    // With fixed ions, each classification covers the active part of its molecule, whose atoms
    // map back to the molecule's atoms through the partition.
    const auto* partition = AssessmentAccess::partition(assessment);
    if (classifications.size() != molecules.size()) {
        throw std::logic_error{"the plan does not classify every molecule of its assessment"};
    }
    auto result = std::vector<std::vector<std::optional<std::size_t>>>{};
    result.reserve(molecules.size());
    for (std::size_t index = 0; index < molecules.size(); ++index) {
        auto& entries = result.emplace_back(molecules[index].atom_count());
        const auto atoms = classifications[index].atom().parameter_entry_indices();
        if (atoms.empty()) {
            continue;
        }
        if (partition == nullptr) {
            if (atoms.size() != entries.size()) {
                throw std::logic_error{"an atom classification does not cover its molecule"};
            }
            std::ranges::copy(atoms, entries.begin());
            continue;
        }
        const auto& active = partition->targets.at(index).active_atom_indices;
        if (atoms.size() != active.size()) {
            throw std::logic_error{"an atom classification does not cover the active atoms"};
        }
        for (std::size_t atom = 0; atom < atoms.size(); ++atom) {
            entries.at(active[atom]) = atoms[atom];
        }
    }
    return result;
}

auto single_bond_entries(const calculation::AssessmentResult& assessment,
                         const calculation::ExecutionPlan& plan) -> std::vector<SingleBondEntry> {
    const auto classifications = AssessmentAccess::classifications(assessment, plan);
    const auto& set = parameter_set(plan);
    const auto& molecules = assessment.molecules();
    constexpr auto none = std::numeric_limits<std::size_t>::max();

    auto matched = std::vector<bool>(set.atom().size(), false);
    auto single_bonded = std::vector<bool>(set.atom().size(), true);
    // The atom entry of the one-bond atom of every bond of each bond entry, or none; a bond entry
    // with bonds that disagree follows no single atom entry.
    auto follows = std::vector<std::optional<std::size_t>>(set.bond().size());
    auto mixed = std::vector<bool>(set.bond().size(), false);
    auto bonds_of = std::vector<std::vector<std::size_t>>(set.atom().size());

    for (std::size_t index = 0; index < classifications.size(); ++index) {
        const auto& molecule = molecules.at(index);
        const auto atoms = classifications[index].atom().parameter_entry_indices();
        const auto bonds = classifications[index].bond().parameter_entry_indices();
        // With fixed ions, the classification covers the calculated part, not the molecule.
        const auto aligned = atoms.size() == molecule.atom_count() &&
                             (bonds.empty() || bonds.size() == molecule.bond_count());
        auto degree = std::vector<std::size_t>(molecule.atom_count(), 0);
        for (std::size_t bond = 0; bond < molecule.bond_count(); ++bond) {
            ++degree[molecule.bond(bond).first_atom_index()];
            ++degree[molecule.bond(bond).second_atom_index()];
        }
        for (std::size_t atom = 0; atom < atoms.size(); ++atom) {
            matched[atoms[atom]] = true;
            if (!aligned || degree[atom] != 1) {
                single_bonded[atoms[atom]] = false;
            }
        }
        if (!aligned) {
            continue;
        }
        for (std::size_t bond = 0; bond < molecule.bond_count(); ++bond) {
            const auto first = molecule.bond(bond).first_atom_index();
            const auto second = molecule.bond(bond).second_atom_index();
            const auto ends =
                static_cast<int>(degree[first] == 1) + static_cast<int>(degree[second] == 1);
            if (ends == 2) {
                // Both atoms have one bond, as in HCl, also without bond parameters.
                single_bonded[atoms[first]] = false;
                single_bonded[atoms[second]] = false;
            }
            if (bonds.empty()) {
                continue;
            }
            const auto end = ends != 1 ? none : degree[first] == 1 ? atoms[first] : atoms[second];
            const auto entry = bonds[bond];
            if (end != none) {
                bonds_of[end].push_back(entry);
            }
            if (!follows[entry].has_value()) {
                follows[entry] = end;
            } else if (*follows[entry] != end) {
                mixed[entry] = true;
            }
        }
    }

    auto result = std::vector<SingleBondEntry>(set.atom().size());
    for (std::size_t entry = 0; entry < result.size(); ++entry) {
        auto own = std::move(bonds_of[entry]);
        std::ranges::sort(own);
        own.erase(std::ranges::unique(own).begin(), own.end());
        result[entry].single_bonded = matched[entry] && single_bonded[entry];
        result[entry].exclusive_bonds = std::ranges::all_of(
            own, [&](const std::size_t bond) { return !mixed[bond] && follows[bond] == entry; });
        result[entry].bond_entries = std::move(own);
    }
    return result;
}

} // namespace chargefw::parameterization
