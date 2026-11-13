#pragma once

#include <chargefw/calculation/charge_issues.h>
#include <chargefw/charges/charge_set.h>
#include <chargefw/core/molecule_collection.h>

#include "calculation/fixed_charge_partition.h"

#include <vector>

namespace chargefw::calculation::detail {

// The charge issues of the charges, in assignment order. molecules are the original molecules
// the charges belong to; the charges of fixed ions, which are prescribed, are not checked.
[[nodiscard]] auto find_charge_issues(const core::MoleculeCollection& molecules,
                                      const charges::ChargeSet& charges,
                                      const FixedChargePartition* fixed_charge_partition)
    -> std::vector<ChargeIssue>;

} // namespace chargefw::calculation::detail
