#pragma once

#include "methods/parameter_columns.h"

#include <chargefw/methods/calculation_input.h>
#include <chargefw/methods/method_derivatives.h>

#include <cstddef>
#include <memory>
#include <span>
#include <vector>

namespace chargefw::methods::builtin::sqe_core {

// The gradient of a weighted sum of the charges with respect to the SQE parameters, and with
// respect to the initial charges.
struct Gradients {
    ParameterGradient parameters;
    std::vector<double> initial_charges;
};

// The charges of calculate() with their derivatives (sqe_derivatives.cpp). gradient() has the atom
// channels electronegativity, hardness, and width, and the bond channel kappa; it treats the
// initial charges as constant.
class Linearization : public LinearizedCharges {
  public:
    // gradient(), together with the gradient with respect to the initial charges, for variants
    // whose initial charges depend on parameters.
    [[nodiscard]] virtual auto gradients(std::span<const double> cotangent) const -> Gradients = 0;

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
