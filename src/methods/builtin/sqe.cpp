#include "methods/builtin/sqe.h"
#include "methods/atom_parameter_values.h"
#include "methods/builtin/prerequisite_helpers.h"
#include "methods/builtin/sqe_system.h"

#include "methods/fixed_source_validation.h"

#include <chargefw/core/molecule.h>
#include <chargefw/core/position.h>
#include <chargefw/parameters/models/parameter_view.h>

#include <Eigen/Cholesky>
#include <Eigen/Core>

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace chargefw::methods::builtin {
namespace {

// erf(x) rounds to exactly 1.0 for x >= 6, so pairs with distance >= 6 * sqrt(width_sum) reduce to
// 1/r without evaluating erf. Comparing squares avoids the sqrt and also covers zero widths.
[[nodiscard]] auto interaction(const double distance, const double width_i, const double width_j)
    -> double {
    const auto width_sum = 2.0 * width_i * width_i + 2.0 * width_j * width_j;

    if (distance * distance >= 36.0 * width_sum) {
        return 1.0 / distance;
    }

    return std::erf(distance / std::sqrt(width_sum)) / distance;
}

} // namespace

auto SQEMethod::add_method_specific_prerequisite_issues(const MethodPrerequisiteInput& input,
                                                        PrerequisiteResult& result) const -> void {
    detail::add_component_neutrality_prerequisite_issue(
        input, result, "SQE", "it moves charge only along bonds, starting from zero charges");
}

auto sqe_core::estimate(const TargetSize& target, const bool initial_charges) -> ResourceEstimate {
    const auto n = target.atoms;
    const auto m = target.bonds;
    if (m == 0.0) {
        // calculate() returns zero or initial charges without building either system.
        return {.memory_bytes = 0.0, .operations = n};
    }
    // Dense n x n atom interactions and the m x m split-charge system factorized in place by
    // Cholesky. Initial charges add one dense matrix-vector product.
    return {.memory_bytes = 8.0 * n * n + 8.0 * m * m,
            .operations = 1.0 / 3.0 * m * m * m + 0.5 * n * n + m * m +
                          (initial_charges ? 2.0 * n * n : 0.0)};
}

auto sqe_core::build_system(const CalculationInput& input,
                            const std::span<const double> initial_charge_values) -> System {
    const auto& molecule = input.molecule();
    const auto& geometry = input.geometry();
    const auto& parameters = input.parameters();
    const auto fixed_sources = input.fixed_sources();

    const auto atom_count = molecule.atom_count();
    const auto bond_count = molecule.bond_count();

    const auto electronegativity = parameters.atom("electronegativity");
    const auto hardness = parameters.atom("hardness");
    // The pair loop below reads widths O(n^2) times; resolve them once per atom.
    const auto width = methods::detail::atom_parameter_values(parameters, "width", atom_count);
    const auto kappa = parameters.bond("kappa");

    const auto n = static_cast<Eigen::Index>(atom_count);
    const auto m = static_cast<Eigen::Index>(bond_count);

    auto system = System{};
    auto& charge_matrix = system.charge_matrix;
    auto& initial_charges = system.initial_charges;
    charge_matrix = Eigen::MatrixXd::Zero(n, n);
    Eigen::VectorXd charge_rhs = Eigen::VectorXd::Zero(n);

    if (!initial_charge_values.empty()) {
        initial_charges.resize(n);
        for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
            initial_charges(static_cast<Eigen::Index>(atom_index)) =
                initial_charge_values[atom_index];
        }
    }

    for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
        const auto i = static_cast<Eigen::Index>(atom_index);
        charge_matrix(i, i) = hardness[atom_index];
        charge_rhs(i) = -electronegativity[atom_index];
        const auto& position = geometry.position(atom_index);
        for (const auto& source : fixed_sources) {
            const auto distance = core::distance(position, source.position);
            charge_rhs(i) -= source.charge * interaction(distance, width[atom_index], 0.0);
        }

        for (std::size_t other_atom_index = atom_index + 1; other_atom_index < atom_count;
             ++other_atom_index) {
            const auto j = static_cast<Eigen::Index>(other_atom_index);
            const auto value = interaction(geometry.distance(atom_index, other_atom_index),
                                           width[atom_index], width[other_atom_index]);
            charge_matrix(i, j) = value;
            charge_matrix(j, i) = value;
        }
    }

    if (!initial_charge_values.empty()) {
        charge_rhs -= charge_matrix * initial_charges;
        charge_rhs += charge_matrix.diagonal().cwiseProduct(initial_charges);
    }

    for (std::size_t bond_index = 0; bond_index < bond_count; ++bond_index) {
        const auto& bond = molecule.bond(bond_index);
        system.bonds.emplace_back(static_cast<Eigen::Index>(bond.first_atom_index()),
                                  static_cast<Eigen::Index>(bond.second_atom_index()));
    }

    // Split-charge system S p = T r with S = T J T^T + diag(kappa), where the transfer matrix T has
    // row a = +1 at atom i and -1 at atom j for bond a = (i, j). T is not formed; for bonds
    // a = (i, j) and b = (k, l):
    //   S(a, b) = J(i, k) - J(i, l) - J(j, k) + J(j, l)  (+ kappa(a) if a == b)
    //   (T r)(a) = r(i) - r(j)
    // Atom charges are q = T^T p: q(i) += p(a), q(j) -= p(a).
    auto& split_matrix = system.split_matrix;
    auto& split_rhs = system.split_rhs;
    split_matrix.resize(m, m);
    split_rhs.resize(m);
    for (Eigen::Index a = 0; a < m; ++a) {
        const auto [i, j] = system.bonds[static_cast<std::size_t>(a)];
        split_rhs(a) = charge_rhs(i) - charge_rhs(j);
        for (Eigen::Index b = a; b < m; ++b) {
            const auto [k, l] = system.bonds[static_cast<std::size_t>(b)];
            const auto value = charge_matrix(i, k) - charge_matrix(i, l) - charge_matrix(j, k) +
                               charge_matrix(j, l);
            split_matrix(a, b) = value;
            split_matrix(b, a) = value;
        }
        split_matrix(a, a) += kappa[static_cast<std::size_t>(a)];
    }

    return system;
}

auto sqe_core::solve(System& system) -> void {
    // S is symmetric; a positive definite S means the energy has a unique minimum. Any other
    // system (indefinite, semidefinite, or nearly singular to working precision) is rejected
    // rather than solved. The calling execution adds the molecule, conformer, and fragment.
    const Eigen::LLT<Eigen::Ref<Eigen::MatrixXd>> llt{system.split_matrix};
    if (llt.info() != Eigen::Success) {
        throw std::runtime_error{
            "SQE split-charge system could not be factorized as positive definite; no charges were "
            "returned. Check the geometry or use another parameter set"};
    }
    system.split_charges = llt.solve(system.split_rhs);
}

auto sqe_core::atom_charges(const System& system) -> Eigen::VectorXd {
    Eigen::VectorXd charges = Eigen::VectorXd::Zero(system.charge_matrix.rows());
    for (Eigen::Index a = 0; a < system.split_charges.size(); ++a) {
        const auto [i, j] = system.bonds[static_cast<std::size_t>(a)];
        charges(i) += system.split_charges(a);
        charges(j) -= system.split_charges(a);
    }

    if (system.initial_charges.size() != 0) {
        charges += system.initial_charges;
    }

    return charges;
}

auto sqe_core::validate_input(const CalculationInput& input,
                              const std::span<const double> initial_charge_values) -> void {
    methods::detail::validate_fixed_source_input(input, "SQE-family");
    const auto atom_count = input.molecule().atom_count();
    if (!initial_charge_values.empty() && initial_charge_values.size() != atom_count) {
        throw std::logic_error{"SQE-family initial charge count must match atom count"};
    }
}

auto sqe_core::calculate(const CalculationInput& input,
                         const std::span<const double> initial_charge_values)
    -> std::vector<double> {
    validate_input(input, initial_charge_values);
    const auto& molecule = input.molecule();

    const auto atom_count = molecule.atom_count();
    const auto bond_count = molecule.bond_count();

    if (atom_count == 0) {
        return {};
    }

    if (bond_count == 0) {
        if (initial_charge_values.empty()) {
            return std::vector<double>(atom_count, 0.0);
        }

        return {initial_charge_values.begin(), initial_charge_values.end()};
    }

    auto system = build_system(input, initial_charge_values);
    solve(system);
    const Eigen::VectorXd charges = atom_charges(system);
    return {charges.data(), charges.data() + charges.size()};
}

auto SQEMethod::calculate(const CalculationInput& input) const -> charges::AtomicCharges {
    return charges::AtomicCharges{sqe_core::calculate(input)};
}

auto SQEMethod::estimate(const TargetSize& target, const MethodOptions& /*options*/) const
    -> ResourceEstimate {
    return sqe_core::estimate(target, false);
}

} // namespace chargefw::methods::builtin
