#include "methods/builtin/sqeqp.h"

#include "methods/builtin/sqe.h"
#include "methods/reference_charges.h"

#include <chargefw/core/molecule.h>
#include <chargefw/parameters/models/parameter_view.h>

#include <cstddef>
#include <utility>
#include <vector>

namespace chargefw::methods::builtin {

auto SQEqpMethod::calculate(const CalculationInput& input) const -> charges::AtomicCharges {
    const auto& molecule = input.molecule();
    const auto& parameters = input.parameters();

    const auto atom_count = molecule.atom_count();

    if (atom_count == 0) {
        return charges::AtomicCharges{sqe_core::calculate(input)};
    }

    if (const auto references = input.reference_charges(); !references.empty()) {
        return charges::AtomicCharges{sqe_core::calculate(input, references)};
    }

    const auto initial_charge_parameter = parameters.atom("q0");
    auto initial_charges = std::vector<double>(atom_count);
    for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
        initial_charges[atom_index] = initial_charge_parameter[atom_index];
    }
    initial_charges = methods::detail::component_normalized_charges(
        molecule, input.topology(), std::move(initial_charges), input.target_charge());

    return charges::AtomicCharges{sqe_core::calculate(input, initial_charges)};
}

auto SQEqpMethod::estimate(const TargetSize& target, const MethodOptions& /*options*/) const
    -> ResourceEstimate {
    return sqe_core::estimate(target, true);
}

} // namespace chargefw::methods::builtin
