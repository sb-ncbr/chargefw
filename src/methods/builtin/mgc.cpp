#include "methods/builtin/mgc.h"
#include "methods/linear_algebra/dense_solve.h"

#include "methods/builtin/prerequisite_helpers.h"

#include <chargefw/core/periodic_table.h>

#include <Eigen/Core>

#include <cmath>
#include <cstddef>
#include <vector>

namespace chargefw::methods::builtin {
auto MGCMethod::add_method_specific_prerequisite_issues(const MethodPrerequisiteInput& input,
                                                        PrerequisiteResult& result) const -> void {
    detail::add_element_prerequisite_issues(
        input, result, "MGC requires a positive element electronegativity",
        [](const core::Element& element) -> bool { return element.electronegativity > 0.0; });

    detail::add_component_neutrality_prerequisite_issue(input, result, "MGC");
}

auto MGCMethod::calculate(const CalculationInput& input) const -> charges::AtomicCharges {
    const auto& molecule = input.molecule();
    const auto& table = core::periodic_table();

    const auto n = static_cast<Eigen::Index>(molecule.atom_count());

    if (n == 0) {
        return charges::AtomicCharges{std::vector<double>{}};
    }

    Eigen::MatrixXd S = Eigen::MatrixXd::Zero(n, n);
    Eigen::VectorXd X0 = Eigen::VectorXd::Zero(n);

    auto log_sum = 0.0;

    for (Eigen::Index i = 0; i < n; ++i) {
        const auto& atom = molecule.atom(static_cast<std::size_t>(i));
        const auto& element = table.element(atom.atomic_number());

        X0(i) = element.electronegativity;
        S(i, i) = 1.0;
        log_sum += std::log(X0(i));
    }

    for (const auto& bond : molecule.bonds()) {
        const auto i = static_cast<Eigen::Index>(bond.first_atom_index());
        const auto j = static_cast<Eigen::Index>(bond.second_atom_index());
        const auto w = static_cast<double>(core::bond_order_value(bond.order()));

        S(i, i) += w;
        S(j, j) += w;
        S(i, j) -= w;
        S(j, i) -= w;
    }

    Eigen::VectorXd chi = linear_algebra::solve_in_place(S, X0);
    chi -= X0;
    chi /= std::exp(log_sum / static_cast<double>(n));

    return charges::AtomicCharges{std::vector<double>{chi.data(), chi.data() + chi.size()}};
}

auto MGCMethod::estimate(const TargetSize& target, const MethodOptions& /*options*/) const
    -> ResourceEstimate {
    const auto n = target.atoms;
    // Dense n x n system factorized by LU.
    return {.memory_bytes = 8.0 * n * n, .operations = 2.0 / 3.0 * n * n * n};
}

} // namespace chargefw::methods::builtin
