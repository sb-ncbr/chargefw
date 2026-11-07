#pragma once

#include <chargefw/methods/calculation_input.h>

#include <Eigen/Core>

namespace chargefw::methods::builtin::eem_core {

// The EEM equations of one input with at least one atom, M x = b with x = (q, chi_equalized):
//
//   M = [ A  1 ]   A_ii = B_i,   A_ij = kappa / R_ij,
//       [ 1  0 ]   b_i = -A_i - kappa sum_s Q_s / r_is,   b_n = total charge.
//
// Shared by calculate() and linearize(), so both solve the same system.
struct System {
    Eigen::MatrixXd matrix;
    Eigen::VectorXd rhs;
};

[[nodiscard]] auto build_system(const CalculationInput& input) -> System;

} // namespace chargefw::methods::builtin::eem_core
