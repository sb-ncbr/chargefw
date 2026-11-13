#pragma once

#include "methods/method_applicability.h"

#include <chargefw/calculation/execution_policy.h>
#include <chargefw/calculation/observer.h>
#include <chargefw/charges/charge_collection.h>
#include <chargefw/features/prepared_molecule_collection.h>

#include <cstddef>

namespace chargefw::calculation {
namespace detail {
struct FixedChargePartition;

// Radius within which a cover pivot retains atoms during cover execution.
inline constexpr double cover_retained_radius = 3.0;
} // namespace detail

[[nodiscard]] auto calculate_cover_charges(
    const methods::ApplicableMethod& selected,
    const features::PreparedMoleculeCollection& molecules, const ExecutionPolicy& policy,
    std::size_t max_threads, const CalculationObserver& observer,
    const detail::FixedChargePartition* fixed_charge_partition = nullptr) -> charges::ChargeSet;

} // namespace chargefw::calculation
