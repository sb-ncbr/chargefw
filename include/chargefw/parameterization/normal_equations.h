#pragma once

#include <chargefw/calculation/assessment.h>
#include <chargefw/methods/method.h>
#include <chargefw/parameters/models/parameter_set.h>

#include <cstddef>
#include <span>
#include <vector>

namespace chargefw::parameterization {

// The quantities of a weighted least-squares fit of charges at one set of parameter values. Rows
// are every charge of every assignment, concatenated in assignment order, fixed ions included;
// slots are those of parameter_layout() (parameter_vector.h) of the candidate. With residuals
// r = charges - references, weights w, and the Jacobian J of the charges with respect to the slots:
struct NormalEquations {
    // One per row.
    std::vector<double> charges;
    // sum_i (w_i r_i)^2.
    double sum_of_squares = 0.0;
    // J^T diag(w^2) J, slots x slots, row-major: the Gauss-Newton matrix, not the exact Hessian.
    std::vector<double> gauss_newton;
    // J^T diag(w^2) r, one per slot.
    std::vector<double> gradient;
};

// Whether normal_equations() has exact derivatives for `method`: EEM, SQE, SQE+q0, and SQE+qp.
[[nodiscard]] auto has_exact_derivatives(const methods::Method& method) -> bool;

// The normal equations of `plan` with the values of `candidate`, from the exact derivatives of the
// plan's method. `references` and `weights` have one finite value per row; weights are not
// negative, and a row with weight zero is not used. Rows of fixed ions have zero derivatives but
// count in the sum of squares. For identical inputs in the same executable and numerical
// environment, the result does not depend on `max_threads`. The calculation cannot be cancelled.
//
// Throws std::invalid_argument for a candidate that calculation::calculate(plan, candidate) would
// reject, a plan without full execution, a method without exact derivatives, or malformed
// references or weights. A failing calculation, or one whose contribution is not finite, throws
// std::runtime_error or std::logic_error naming the molecule and conformer: the failures
// ExecutionResult reports as numerical_failure. Finite contributions whose sum is not finite throw
// std::runtime_error without naming a molecule.
[[nodiscard]] auto normal_equations(const calculation::ExecutionPlan& plan,
                                    const parameters::ParameterSet& candidate,
                                    std::span<const double> references,
                                    std::span<const double> weights, std::size_t max_threads = 0)
    -> NormalEquations;

} // namespace chargefw::parameterization
