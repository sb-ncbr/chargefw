#pragma once

#include <chargefw/methods/method_requirements.h>

#include <Eigen/Core>
#include <Eigen/LU>

namespace chargefw::methods::builtin {

// The final row and column are reserved for the total-charge constraint.
inline auto set_charge_constraint(Eigen::MatrixXd& matrix, Eigen::VectorXd& rhs,
                                  const double target_charge) -> void {
    const auto n = matrix.rows() - 1;
    matrix.row(n).setOnes();
    matrix.col(n).setOnes();
    matrix(n, n) = 0.0;
    rhs(n) = target_charge;
}

// Solves matrix * x = rhs with partial-pivoting LU. The matrix is overwritten by its factors
// instead of copied, which halves peak memory and cache footprint for large dense systems.
[[nodiscard]] inline auto solve_in_place(Eigen::MatrixXd& matrix, const Eigen::VectorXd& rhs)
    -> Eigen::VectorXd {
    const Eigen::PartialPivLU<Eigen::Ref<Eigen::MatrixXd>> lu{matrix};
    return lu.solve(rhs);
}

// Returns the full solution; the last entry is the Lagrange multiplier. Callers take the charges
// with head(n) or map them as ABEEM does.
[[nodiscard]] inline auto solve_charge_constrained(Eigen::MatrixXd& matrix, Eigen::VectorXd& rhs,
                                                   const double target_charge) -> Eigen::VectorXd {
    set_charge_constraint(matrix, rhs, target_charge);
    return solve_in_place(matrix, rhs);
}

[[nodiscard]] inline auto charge_constrained_estimate(const double n) -> ResourceEstimate {
    const auto size = n + 1.0;
    // Dense (n + 1) x (n + 1) solve plus all-pairs matrix assembly.
    return {.memory_bytes = 8.0 * size * size,
            .operations = 2.0 / 3.0 * size * size * size + 0.5 * n * n};
}

} // namespace chargefw::methods::builtin
