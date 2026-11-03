#pragma once

#include <chargefw/calculation/input_issues.h>
#include <chargefw/core/molecule_collection.h>

#include <vector>

namespace chargefw::calculation::detail {

// The input issues of the molecules, in molecule order.
[[nodiscard]] auto find_input_issues(const core::MoleculeCollection& molecules)
    -> std::vector<InputIssue>;

} // namespace chargefw::calculation::detail
