#pragma once

#include <chargefw/methods/calculation_input.h>

#include <Eigen/Core>

#include <span>
#include <utility>
#include <vector>

namespace chargefw::methods::builtin::sqe_core {

// The SQE equations of one input with at least one bond. Shared by calculate() and by the exact
// derivatives (src/parameterisation), so that both solve the same system.
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

// The initial charges of SQE+q0 (sqeq0.cpp): the formal charges.
[[nodiscard]] auto formal_initial_charges(const CalculationInput& input) -> std::vector<double>;

// The initial charges of SQE+qp (sqeqp.cpp), as its calculate() uses them.
struct ParameterizedInitialCharges {
    // None without atoms; the reference charges when the input has them; otherwise the atom
    // parameter q0, normalized per bonded part.
    std::vector<double> values;
    // Whether `values` come from q0, rather than from reference charges or no atoms.
    bool from_q0 = false;
};
[[nodiscard]] auto parameterized_initial_charges(const CalculationInput& input)
    -> ParameterizedInitialCharges;

} // namespace chargefw::methods::builtin::sqe_core
