#include "methods/builtin/sqeqp.h"

#include "methods/builtin/sqe.h"
#include "methods/builtin/sqe_system.h"
#include "methods/reference_charges.h"

#include <chargefw/core/molecule.h>
#include <chargefw/parameters/models/parameter_view.h>

#include <cstddef>
#include <utility>
#include <vector>

namespace chargefw::methods::builtin {

auto sqe_core::parameterized_initial_charges(const CalculationInput& input)
    -> ParameterizedInitialCharges {
    const auto& molecule = input.molecule();
    const auto atom_count = molecule.atom_count();

    if (atom_count == 0) {
        return {};
    }

    if (const auto references = input.reference_charges(); !references.empty()) {
        return {.values = {references.begin(), references.end()}, .from_q0 = false};
    }

    const auto initial_charge_parameter = input.parameters().atom("q0");
    auto initial_charges = std::vector<double>(atom_count);
    for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
        initial_charges[atom_index] = initial_charge_parameter[atom_index];
    }
    return {.values = methods::detail::component_normalized_charges(
                molecule, input.topology(), std::move(initial_charges), input.target_charge()),
            .from_q0 = true};
}

auto SQEqpMethod::calculate(const CalculationInput& input) const -> charges::AtomicCharges {
    return charges::AtomicCharges{
        sqe_core::calculate(input, sqe_core::parameterized_initial_charges(input).values)};
}

auto SQEqpMethod::estimate(const TargetSize& target, const MethodOptions& /*options*/) const
    -> ResourceEstimate {
    return sqe_core::estimate(target, true);
}

} // namespace chargefw::methods::builtin
