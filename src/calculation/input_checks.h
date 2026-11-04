#pragma once

#include <chargefw/calculation/input_issues.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/features/prepared_molecule_collection.h>

#include "calculation/fixed_charge_partition.h"

#include <span>

#include <vector>

namespace chargefw::calculation::detail {

// The input issues of the molecules, in molecule order.
[[nodiscard]] auto find_input_issues(const core::MoleculeCollection& molecules)
    -> std::vector<InputIssue>;

// The assessment warnings, which are in molecule order, with an atoms_too_close issue added for
// each conformer of the molecules that has close atoms. prepared holds the molecules, or with fixed
// ions their active parts; their cached geometry defects are read. A molecule with fixed ions is
// searched anew instead, so that the ions are included and atom indices are original.
[[nodiscard]] auto add_close_atom_issues(std::span<const InputIssue> assessment_warnings,
                                         const core::MoleculeCollection& molecules,
                                         const features::PreparedMoleculeCollection& prepared,
                                         const FixedChargePartition* fixed_charge_partition)
    -> std::vector<InputIssue>;

} // namespace chargefw::calculation::detail
