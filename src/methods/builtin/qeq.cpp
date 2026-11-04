#include "methods/builtin/qeq.h"
#include "methods/builtin/dense_solve.h"

#include <chargefw/core/molecule.h>
#include <chargefw/parameters/models/parameter_view.h>

#include <Eigen/Core>

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace chargefw::methods::builtin {
namespace {

constexpr auto coulomb = 14.4;
constexpr auto bohr_radius_angstrom = 0.529177;

using OverlapTerm = double (*)(double distance, double hardness_i, double hardness_j);

[[nodiscard]] auto resolve_overlap_term(const std::string& type) -> OverlapTerm {
    if (type == "Nishimoto-Mataga") {
        return [](const double distance, const double hardness_i, const double hardness_j) {
            return coulomb / (distance + 2.0 * coulomb / (hardness_i + hardness_j));
        };
    }
    if (type == "Nishimoto-Mataga-Weiss") {
        return [](const double distance, const double hardness_i, const double hardness_j) {
            constexpr auto factor = 1.2;
            return factor * coulomb /
                   (distance + 2.0 * factor * coulomb / (hardness_i + hardness_j));
        };
    }
    if (type == "Ohno") {
        return [](const double distance, const double hardness_i, const double hardness_j) {
            return coulomb / std::sqrt(distance * distance +
                                       std::pow(2.0 * coulomb / (hardness_i + hardness_j), 2.0));
        };
    }
    if (type == "Ohno-Klopman") {
        return [](const double distance, const double hardness_i, const double hardness_j) {
            return coulomb /
                   std::sqrt(
                       distance * distance +
                       std::pow(coulomb / (2.0 * hardness_i) + coulomb / (2.0 * hardness_j), 2.0));
        };
    }
    if (type == "DasGupta-Huzinaga") {
        return [](const double distance, const double hardness_i, const double hardness_j) {
            constexpr auto k = 0.4 / bohr_radius_angstrom;
            return coulomb / (distance + coulomb / (hardness_i / 2.0 * std::exp(k * distance) +
                                                    hardness_j / 2.0 * std::exp(k * distance)));
        };
    }
    throw std::logic_error{"Unsupported QEq overlap term '" + type + "'"};
}

} // namespace

auto QEqMethod::calculate(const CalculationInput& input) const -> charges::AtomicCharges {
    const auto& molecule = input.molecule();
    const auto& parameters = input.parameters();
    const auto& geometry = input.geometry();

    const auto atom_count = molecule.atom_count();

    if (atom_count == 0) {
        return charges::AtomicCharges{std::vector<double>{}};
    }

    const auto electronegativity = parameters.atom("electronegativity");
    // The pair loop below reads hardness O(n^2) times; resolve it once per atom.
    const auto hardness_parameter = parameters.atom("hardness");
    std::vector<double> hardness(atom_count);
    for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
        hardness[atom_index] = hardness_parameter[atom_index];
    }
    const auto& type = input.method_options().get<std::string>("overlap_term");

    for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
        if (hardness[atom_index] <= 0.0) {
            throw std::logic_error{"QEq atom parameter 'hardness' must be positive"};
        }
    }

    const auto overlap_term = resolve_overlap_term(type);

    const auto n = static_cast<Eigen::Index>(atom_count);

    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n + 1, n + 1);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(n + 1);

    for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
        const auto i = static_cast<Eigen::Index>(atom_index);

        A(i, i) = hardness[atom_index];
        b(i) = -electronegativity[atom_index];

        for (std::size_t other_atom_index = atom_index + 1; other_atom_index < atom_count;
             ++other_atom_index) {
            const auto j = static_cast<Eigen::Index>(other_atom_index);
            const auto value = overlap_term(geometry.distance(atom_index, other_atom_index),
                                            hardness[atom_index], hardness[other_atom_index]);

            A(i, j) = value;
            A(j, i) = value;
        }
    }

    A.row(n).setOnes();
    A.col(n).setOnes();
    A(n, n) = 0.0;
    b(n) = input.target_charge();

    const Eigen::VectorXd q = solve_in_place(A, b).head(n);

    return charges::AtomicCharges{std::vector<double>{q.data(), q.data() + q.size()}};
}

auto QEqMethod::estimate(const TargetSize& target, const MethodOptions& /*options*/) const
    -> ResourceEstimate {
    const auto n = target.atoms;
    const auto size = n + 1.0;
    // Dense (n + 1) x (n + 1) system with a charge constraint.
    return {.memory_bytes = 8.0 * size * size,
            .operations = 2.0 / 3.0 * size * size * size + 0.5 * n * n};
}

} // namespace chargefw::methods::builtin
