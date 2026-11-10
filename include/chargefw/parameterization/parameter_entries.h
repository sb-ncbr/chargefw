#pragma once

#include <chargefw/calculation/assessment.h>

#include <cstddef>
#include <optional>
#include <vector>

namespace chargefw::parameterization {

// How many calculated atoms each atom entry, and bonds each bond entry, of a plan's parameter set
// matches, as classified during assessment. With fixed ions, their atoms are not counted.
struct EntryCounts {
    std::vector<std::size_t> atoms; // [atom entry]
    std::vector<std::size_t> bonds; // [bond entry]
};

// Throws std::invalid_argument for a plan of another assessment or a method without parameters.
[[nodiscard]] auto entry_counts(const calculation::AssessmentResult& assessment,
                                const calculation::ExecutionPlan& plan) -> EntryCounts;

// The atom entry matched to every atom of every molecule of the assessment, in the molecules' own
// atom order; std::nullopt for an atom that no entry matches: a fixed ion, or every atom when the
// parameter set has no atom entries. Throws as entry_counts().
[[nodiscard]] auto atom_entries(const calculation::AssessmentResult& assessment,
                                const calculation::ExecutionPlan& plan)
    -> std::vector<std::vector<std::optional<std::size_t>>>;

// What the molecules show about the atoms of one atom entry and their bonds.
struct SingleBondEntry {
    // The entry matches at least one atom, every atom it matches has exactly one bond, and none of
    // them is bonded to another atom with one bond (as in HCl). False when the entry matches an
    // atom of a molecule calculated with fixed ions, whose atoms do not line up with the
    // classification.
    bool single_bonded = false;
    // The bond entries of the bonds at those atoms, ascending.
    std::vector<std::size_t> bond_entries;
    // Every bond each of `bond_entries` matches, in molecules calculated without fixed ions, has
    // exactly one atom with one bond, and that atom matches this entry.
    bool exclusive_bonds = false;
};

// One per atom entry of the plan's parameter set. Such an atom exchanges charge only through its
// bond; whether that lets a fit hold one of its parameters is the fitting profile's decision.
// Throws as entry_counts().
[[nodiscard]] auto single_bond_entries(const calculation::AssessmentResult& assessment,
                                       const calculation::ExecutionPlan& plan)
    -> std::vector<SingleBondEntry>;

} // namespace chargefw::parameterization
