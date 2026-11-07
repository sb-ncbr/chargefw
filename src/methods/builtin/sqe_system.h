#pragma once

#include <chargefw/methods/calculation_input.h>

#include <Eigen/Core>

#include <span>
#include <utility>
#include <vector>

namespace chargefw::methods::builtin::sqe_core {

// The SQE equations of one input with at least one bond. Shared by calculate() and linearize(), so
// both solve the same system.
struct System {
    // J: hardness on the diagonal, the atom-atom interactions off it.
    Eigen::MatrixXd charge_matrix;
    // Empty without initial charges.
    Eigen::VectorXd initial_charges;
    // The atoms (i, j) of each bond.
    std::vector<std::pair<Eigen::Index, Eigen::Index>> bonds;
    // S = T J T^T + diag(kappa); after solve(), its lower triangle holds the Cholesky factor.
    Eigen::MatrixXd split_matrix;
    // T r.
    Eigen::VectorXd split_rhs;
    // p, set by solve().
    Eigen::VectorXd split_charges;
};

[[nodiscard]] auto build_system(const CalculationInput& input,
                                std::span<const double> initial_charge_values) -> System;

// Factorizes S in place and solves S p = T r. Throws when S is not positive definite.
auto solve(System& system) -> void;

// q = T^T p + q_initial.
[[nodiscard]] auto atom_charges(const System& system) -> Eigen::VectorXd;

} // namespace chargefw::methods::builtin::sqe_core
