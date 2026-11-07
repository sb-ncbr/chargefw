#include "methods/builtin/sqeq0.h"

#include "methods/builtin/sqe.h"
#include "methods/builtin/sqe_derivatives.h"

#include <chargefw/core/molecule.h>

#include <cstddef>
#include <memory>
#include <vector>

namespace chargefw::methods::builtin {

auto SQEq0Method::calculate(const CalculationInput& input) const -> charges::AtomicCharges {
    const auto& molecule = input.molecule();
    const auto atom_count = molecule.atom_count();
    auto initial_charges = std::vector<double>(atom_count, 0.0);

    for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
        initial_charges[atom_index] = molecule.atom(atom_index).formal_charge();
    }

    return charges::AtomicCharges{sqe_core::calculate(input, initial_charges)};
}

auto SQEq0Method::linearize(const CalculationInput& input) const
    -> std::unique_ptr<LinearizedCharges> {
    const auto& molecule = input.molecule();
    auto initial_charges = std::vector<double>(molecule.atom_count(), 0.0);
    for (std::size_t atom_index = 0; atom_index < molecule.atom_count(); ++atom_index) {
        initial_charges[atom_index] = molecule.atom(atom_index).formal_charge();
    }
    // The formal charges are input, not parameters: their derivative is zero.
    return sqe_core::linearize(input, initial_charges);
}

auto SQEq0Method::estimate(const TargetSize& target, const MethodOptions& /*options*/) const
    -> ResourceEstimate {
    return sqe_core::estimate(target, true);
}

} // namespace chargefw::methods::builtin
