#include "methods/builtin/eqeq.h"
#include "methods/builtin/charge_constraint.h"

#include <chargefw/core/molecule.h>
#include <chargefw/core/periodic_table.h>

#include <Eigen/Core>

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace chargefw::methods::builtin {
namespace {

constexpr auto lambda = 1.2;
constexpr auto k = 14.4;
constexpr auto hydrogen_atomic_number = 1;
constexpr auto hydrogen_electron_affinity = -2.0;

[[nodiscard]] auto electron_affinity_for_eqeq(const core::Element& element) -> double {
    if (element.atomic_number == hydrogen_atomic_number) {
        return hydrogen_electron_affinity;
    }

    return element.electron_affinity;
}

[[nodiscard]] auto interaction(const double distance, const double hardness_i,
                               const double hardness_j) -> double {
    const auto a = std::sqrt(hardness_i * hardness_j) / k;
    const auto overlap =
        std::exp(-a * a * distance * distance) * (2.0 * a - a * a * distance - 1.0 / distance);

    return lambda * k / 2.0 * (1.0 / distance + overlap);
}

} // namespace

auto eqeq_core::calculate(const CalculationInput& input, const std::string_view method_name)
    -> std::vector<double> {
    const auto& molecule = input.molecule();
    const auto& geometry = input.geometry();

    const auto atom_count = molecule.atom_count();

    if (atom_count == 0) {
        return {};
    }

    const auto& table = core::periodic_table();
    const auto n = static_cast<Eigen::Index>(atom_count);

    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(n + 1, n + 1);
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(n + 1);
    Eigen::VectorXd hardness = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd electronegativity = Eigen::VectorXd::Zero(n);

    for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
        const auto i = static_cast<Eigen::Index>(atom_index);
        const auto& atom = molecule.atom(atom_index);
        const auto& element = table.element(atom.atomic_number());
        const auto electron_affinity = electron_affinity_for_eqeq(element);

        electronegativity(i) = (element.first_ionization_potential + electron_affinity) / 2.0;
        hardness(i) = element.first_ionization_potential - electron_affinity;
    }

    if ((hardness.array() <= 0.0).any()) {
        throw std::logic_error{std::string{method_name} + " hardness must be positive"};
    }

    for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
        const auto i = static_cast<Eigen::Index>(atom_index);

        matrix(i, i) = hardness(i);
        rhs(i) = -electronegativity(i);

        for (std::size_t other_atom_index = atom_index + 1; other_atom_index < atom_count;
             ++other_atom_index) {
            const auto j = static_cast<Eigen::Index>(other_atom_index);
            const auto distance = geometry.distance(atom_index, other_atom_index);
            const auto value = interaction(distance, hardness(i), hardness(j));

            matrix(i, j) = value;
            matrix(j, i) = value;
        }
    }

    const Eigen::VectorXd q = solve_charge_constrained(matrix, rhs, input.target_charge()).head(n);

    return {q.data(), q.data() + q.size()};
}

auto EQeqMethod::calculate(const CalculationInput& input) const -> charges::AtomicCharges {
    return charges::AtomicCharges{eqeq_core::calculate(input, "EQeq")};
}

auto EQeqMethod::estimate(const TargetSize& target, const MethodOptions& /*options*/) const
    -> ResourceEstimate {
    return charge_constrained_estimate(target.atoms);
}

} // namespace chargefw::methods::builtin
