#pragma once

#include <chargefw/adapters/charge_result.h>
#include <chargefw/charges/charge_collection.h>

#include <cstddef>
#include <span>

namespace chargefw::adapters {

// The charge assignments of one input record: one for the molecule, or one per conformer in
// conformer order. Results keep assignments in input order, so writers find a record's assignments
// by binary search instead of scanning those of every record.
[[nodiscard]] auto record_assignments(const ChargeCalculationResult& result,
                                      std::size_t record_index)
    -> std::span<const charges::ChargeAssignment>;

} // namespace chargefw::adapters
