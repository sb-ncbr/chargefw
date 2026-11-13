#pragma once

#include "methods/method_applicability.h"

#include <chargefw/calculation/observer.h>
#include <chargefw/charges/charge_collection.h>
#include <chargefw/features/prepared_molecule_collection.h>

#include <cstddef>

namespace chargefw::calculation {
namespace detail {
struct FixedChargePartition;
}

[[nodiscard]] auto
calculate_full_charges(const methods::ApplicableMethod& selected,
                       const features::PreparedMoleculeCollection& molecules,
                       std::size_t max_threads, const CalculationObserver& observer,
                       const detail::FixedChargePartition* fixed_charge_partition = nullptr)
    -> charges::ChargeSet;

} // namespace chargefw::calculation
