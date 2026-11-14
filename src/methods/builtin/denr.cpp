#include "methods/builtin/denr.h"

#include "methods/builtin/prerequisite_helpers.h"
#include "methods/linear_algebra/dense_solve.h"

#include <chargefw/parameters/models/parameter_view.h>

#include <stdexcept>
#include <vector>

namespace chargefw::methods::builtin {

auto DENRMethod::add_method_specific_prerequisite_issues(const MethodPrerequisiteInput& input,
                                                         PrerequisiteResult& result) const -> void {
    detail::add_component_neutrality_prerequisite_issue(input, result, "DENR");
}

auto DENRMethod::calculate(const CalculationInput& input) const -> charges::AtomicCharges {
    const auto& molecule = input.molecule();
    const auto& parameters = input.parameters();
    const auto& options = input.method_options();

    const auto n = static_cast<Eigen::Index>(molecule.atom_count());

    if (n == 0) {
        return charges::AtomicCharges{std::vector<double>{}};
    }

    const auto electronegativity = parameters.atom("electronegativity");
    const auto hardness = parameters.atom("hardness");
    const auto dt = options.get<double>("step");
    const auto iterations = options.get<int>("iters");

    // DENR iterates q <- A^-1 (q - dt L chi) with A = I + dt L eta, where L is the bond-graph
    // Laplacian and eta = diag(hardness). Neither L nor eta is formed; for bond (i, j):
    //   A(i, i) += dt eta(i),  A(j, j) += dt eta(j),  A(i, j) -= dt eta(j),  A(j, i) -= dt eta(i)
    //   (dt L chi)(i) += dt (chi(i) - chi(j)),  (dt L chi)(j) -= dt (chi(i) - chi(j))
    Eigen::MatrixXd A = Eigen::MatrixXd::Identity(n, n);
    Eigen::VectorXd tmp = Eigen::VectorXd::Zero(n);

    for (const auto& bond : molecule.bonds()) {
        const auto first = bond.first_atom_index();
        const auto second = bond.second_atom_index();
        const auto i = static_cast<Eigen::Index>(first);
        const auto j = static_cast<Eigen::Index>(second);

        A(i, i) += dt * hardness[first];
        A(j, j) += dt * hardness[second];
        A(i, j) -= dt * hardness[second];
        A(j, i) -= dt * hardness[first];

        const auto transfer = dt * (electronegativity[first] - electronegativity[second]);
        tmp(i) += transfer;
        tmp(j) -= transfer;
    }

    // Factorize in place to avoid holding a second dense copy of A.
    const linear_algebra::InPlaceLu lu{A};
    Eigen::VectorXd q = Eigen::VectorXd::Zero(n);

    for (auto i = 0; i < iterations; ++i) {
        q = lu.solve(Eigen::VectorXd{q - tmp});
    }

    return charges::AtomicCharges{std::vector<double>{q.data(), q.data() + q.size()}};
}

auto DENRMethod::estimate(const TargetSize& target, const MethodOptions& options) const
    -> ResourceEstimate {
    const auto n = target.atoms;
    const auto iterations = options.get<int>("iters");
    // A is dense n x n; LU and one solve per iteration.
    return {.memory_bytes = 8.0 * n * n,
            .operations = 2.0 / 3.0 * n * n * n + iterations * 2.0 * n * n};
}

} // namespace chargefw::methods::builtin
