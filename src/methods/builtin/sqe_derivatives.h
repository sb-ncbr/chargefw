#pragma once

#include "methods/parameter_columns.h"

#include <chargefw/methods/calculation_input.h>
#include <chargefw/methods/method_derivatives.h>

#include <memory>
#include <span>
#include <vector>

namespace chargefw::methods::builtin::sqe_core {

// The charges of calculate() with their derivatives (sqe_derivatives.cpp). jacobian() takes the
// atom parameters electronegativity, hardness, and width, and the bond parameter kappa; it treats
// the initial charges as constant.
class Linearization : public LinearizedCharges {
  public:
    using LinearizedCharges::jacobian;

    // jacobian(), with the initial charges changing too: by column j of initial_charges.values
    // along column initial_charges.columns[j], for variants whose initial charges depend on
    // parameters.
    [[nodiscard]] virtual auto jacobian(const ParameterColumns& columns,
                                        const methods::detail::ColumnBlock& initial_charges) const
        -> std::vector<double> = 0;
};

[[nodiscard]] auto linearize(const CalculationInput& input,
                             std::span<const double> initial_charge_values)
    -> std::unique_ptr<Linearization>;

} // namespace chargefw::methods::builtin::sqe_core
