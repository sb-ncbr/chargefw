#include "methods/builtin/sqeq0.h"

#include "methods/builtin/initial_charges.h"
#include "methods/builtin/sqe.h"

namespace chargefw::methods::builtin {

auto SQEq0Method::calculate(const CalculationInput& input) const -> charges::AtomicCharges {
    return charges::AtomicCharges{
        sqe_core::calculate(input, detail::formal_initial_charges(input.molecule()))};
}

auto SQEq0Method::estimate(const TargetSize& target, const MethodOptions& /*options*/) const
    -> ResourceEstimate {
    return sqe_core::estimate(target, true);
}

} // namespace chargefw::methods::builtin
