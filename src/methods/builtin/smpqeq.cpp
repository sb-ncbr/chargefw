#include "methods/builtin/smpqeq.h"
#include "methods/atom_parameter_values.h"
#include "methods/builtin/charge_constraint.h"

#include <chargefw/core/molecule.h>
#include <chargefw/parameters/models/parameter_view.h>

#include <Eigen/Core>

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace chargefw::methods::builtin {
namespace {

constexpr auto coulomb = 14.4;
constexpr auto maximum_iterations = 100;
constexpr auto convergence_tolerance = 1.0e-10;

[[nodiscard]] auto interaction(const double distance, const double second_i, const double second_j)
    -> double {

    if (second_i <= 0.0 || second_j <= 0.0) {
        throw std::logic_error{"SMP/QEq atom parameter 'second' must be positive"};
    }

    const auto gamma = 2.0 * std::sqrt(second_i * second_j);
    const auto scaled_distance = distance / coulomb;

    return 1.0 / std::cbrt(1.0 / std::pow(gamma, 3.0) + std::pow(scaled_distance, 3.0));
}

[[nodiscard]] auto diagonal_term(const double previous_charge, const double second,
                                 const double third, const double fourth) -> double {
    return 2.0 * (second + third * previous_charge / 4.0 +
                  fourth * previous_charge * previous_charge / 12.0);
}

} // namespace

auto SMPQEqMethod::calculate(const CalculationInput& input) const -> charges::AtomicCharges {
    const auto& molecule = input.molecule();
    const auto& parameters = input.parameters();
    const auto& geometry = input.geometry();

    const auto atom_count = molecule.atom_count();

    if (atom_count == 0) {
        return charges::AtomicCharges{std::vector<double>{}};
    }

    const auto first = parameters.atom("first");
    // The pair loop below reads second O(n^2) times; resolve it once per atom.
    const auto second = methods::detail::atom_parameter_values(parameters, "second", atom_count);
    const auto third = parameters.atom("third");
    const auto fourth = parameters.atom("fourth");

    const auto n = static_cast<Eigen::Index>(atom_count);

    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n + 1, n + 1);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(n + 1);
    Eigen::VectorXd previous_charges = Eigen::VectorXd::Zero(n);

    // Zhang specifies self-consistent iteration, not a fixed iteration count.
    for (auto iteration = 0; iteration < maximum_iterations; ++iteration) {
        A.setZero();

        for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
            const auto i = static_cast<Eigen::Index>(atom_index);
            const auto previous_charge = b(i);

            A(i, i) = diagonal_term(previous_charge, second[atom_index], third[atom_index],
                                    fourth[atom_index]);
            b(i) = -first[atom_index];

            for (std::size_t other_atom_index = atom_index + 1; other_atom_index < atom_count;
                 ++other_atom_index) {
                const auto j = static_cast<Eigen::Index>(other_atom_index);
                const auto distance = geometry.distance(atom_index, other_atom_index);

                const auto value =
                    interaction(distance, second[atom_index], second[other_atom_index]);

                A(i, j) = value;
                A(j, i) = value;
            }
        }

        b = solve_charge_constrained(A, b, input.target_charge());

        if ((b.head(n) - previous_charges).cwiseAbs().maxCoeff() <= convergence_tolerance) {
            return charges::AtomicCharges{std::vector<double>{b.data(), b.data() + n}};
        }

        previous_charges = b.head(n);
    }

    throw std::runtime_error{"SMP/QEq self-consistent solve did not converge"};
}

auto SMPQEqMethod::estimate(const TargetSize& target, const MethodOptions& /*options*/) const
    -> ResourceEstimate {
    auto estimate = charge_constrained_estimate(target.atoms);
    estimate.operations *= maximum_iterations;
    return estimate;
}

} // namespace chargefw::methods::builtin
