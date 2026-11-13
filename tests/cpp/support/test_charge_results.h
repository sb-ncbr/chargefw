#pragma once

#include <chargefw/adapters/charge_result.h>
#include <chargefw/calculation/calculation.h>
#include <chargefw/charges/charge_set.h>

#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace chargefw::test {

/// Effective full-execution provenance for hand-built results.
[[nodiscard]] inline auto full_effective(std::string method_id,
                                         std::optional<std::string> parameter_set_id = {},
                                         std::optional<calculation::FixedIons> fixed_ions = {})
    -> calculation::EffectiveCalculation {
    return {.method_id = std::move(method_id),
            .parameter_set_id = std::move(parameter_set_id),
            .execution_policy = calculation::ExecutionPolicy{},
            .fixed_ions = std::move(fixed_ions)};
}

/// Wraps a charge set as a full-execution result owned by the given records.
[[nodiscard]] inline auto calculation_result(std::vector<adapters::ImportedMoleculeRecord> records,
                                             charges::ChargeSet charge_set)
    -> adapters::ChargeCalculationResult {
    auto method_id = std::string{charge_set.method_id()};
    auto parameter_set_id = charge_set.parameter_set_id().transform(
        [](const std::string_view value) { return std::string{value}; });
    return adapters::make_charge_calculation_result(
        std::move(records), {},
        {.charges = std::move(charge_set),
         .effective = full_effective(std::move(method_id), std::move(parameter_set_id))});
}

} // namespace chargefw::test
