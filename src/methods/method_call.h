#pragma once

#include <chargefw/charges/atomic_charges.h>
#include <chargefw/methods/calculation_input.h>
#include <chargefw/methods/method.h>
#include <chargefw/methods/method_failure.h>

namespace chargefw::methods::detail {

// The failure of a calculation that gave charges that are not finite numbers.
[[nodiscard]] inline auto non_finite_charges_failure() -> MethodFailure {
    return MethodFailure{FailureKind::no_solution,
                         "the calculation gave charges that are not finite numbers, for example "
                         "because its equations have no unique solution; check the geometry or "
                         "use another parameter set"};
}

// Calls the method; charges that are not finite numbers are a no_solution failure.
[[nodiscard]] inline auto calculate_charges(const Method& method, const CalculationInput& input)
    -> charges::AtomicCharges {
    try {
        return method.calculate(input);
    } catch (const charges::NonFiniteCharges&) {
        throw non_finite_charges_failure();
    }
}

} // namespace chargefw::methods::detail
