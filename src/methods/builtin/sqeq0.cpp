#include "methods/builtin/sqeq0.h"

#include "methods/builtin/sqe.h"
#include "methods/builtin/sqe_system.h"

#include <chargefw/core/molecule.h>

#include <cstddef>
#include <vector>

namespace chargefw::methods::builtin {

auto sqe_core::formal_initial_charges(const CalculationInput& input) -> std::vector<double> {
    const auto& molecule = input.molecule();
    const auto atom_count = molecule.atom_count();
    auto initial_charges = std::vector<double>(atom_count, 0.0);

    for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
        initial_charges[atom_index] = molecule.atom(atom_index).formal_charge();
    }

    return initial_charges;
}

auto SQEq0Method::calculate(const CalculationInput& input) const -> charges::AtomicCharges {
    return charges::AtomicCharges{
        sqe_core::calculate(input, sqe_core::formal_initial_charges(input))};
}

auto SQEq0Method::estimate(const TargetSize& target, const MethodOptions& /*options*/) const
    -> ResourceEstimate {
    return sqe_core::estimate(target, true);
}

} // namespace chargefw::methods::builtin
