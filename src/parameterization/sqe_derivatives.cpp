// Exact derivatives of the SQE-family charges.
//
// calculate() (sqe.cpp) solves, for atoms i, j and bonds a:
//
//   S p = T r,   q = T^T p + q_initial,   S = T J T^T + diag(kappa),
//   J_ii = hardness_i,   J_ij = I(d_ij, width_i, width_j),
//   r_i = -electronegativity_i - phi_i(width_i) - sum_{j != i} J_ij q_initial_j,
//
// with phi the fixed sources' potential and the last sum only with initial charges. Differentiating
// both equations, a change of the parameters and initial charges changes the charges by
//
//   S dp = T u - d kappa p,   dq = T^T dp + d q_initial,
//   u = -d electronegativity - d hardness q_split - d width (dphi/dwidth + W q) - W^T (d width q)
//       - J_offdiagonal d q_initial,
//
// with q_split = T^T p, W_ij = dJ_ij/dwidth_i for i != j, and products taken per atom or per bond.
// jacobian() takes one such change per column.

#include "parameterization/method_derivatives.h"
#include "parameterization/parameter_columns.h"

#include "features/topology_helpers.h"
#include "methods/atom_parameter_values.h"
#include "methods/builtin/initial_charges.h"
#include "methods/builtin/sqe.h"
#include "methods/builtin/sqe_system.h"

#include <chargefw/charges/atomic_charges.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/position.h>
#include <chargefw/parameters/models/parameter_view.h>

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <memory>
#include <numbers>
#include <optional>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace chargefw::parameterization::detail {
namespace {

using methods::CalculationInput;
namespace sqe_core = methods::builtin::sqe_core;
namespace builtin_detail = methods::builtin::detail;

// interaction() (sqe.cpp) depends on width_i through s = 2 width_i^2 + 2 width_j^2:
// d/dwidth_i erf(d / sqrt(s)) / d = width_i * g with g = -(4 / sqrt(pi)) exp(-d^2 / s) / s^(3/2).
// Returns g from d^2, which is the same for both widths; it is zero where interaction() returns
// 1/r.
[[nodiscard]] auto interaction_width_factor(const double squared_distance, const double width_i,
                                            const double width_j) -> double {
    const auto width_sum = 2.0 * width_i * width_i + 2.0 * width_j * width_j;

    if (squared_distance >= 36.0 * width_sum) {
        return 0.0;
    }

    return -4.0 * std::numbers::inv_sqrtpi * std::exp(-squared_distance / width_sum) /
           (width_sum * std::sqrt(width_sum));
}

struct WidthDerivatives {
    // dJ_ij/dwidth_i for i != j; zero on the diagonal.
    Eigen::MatrixXd pairs;
    // dphi_i/dwidth_i.
    Eigen::VectorXd sources;
};

[[nodiscard]] auto width_derivatives(const CalculationInput& input) -> WidthDerivatives {
    const auto& geometry = input.geometry();
    const auto atom_count = input.molecule().atom_count();
    const auto n = static_cast<Eigen::Index>(atom_count);
    // As in build_system(): the pair loop reads widths O(n^2) times.
    const auto width =
        methods::detail::atom_parameter_values(input.parameters(), "width", atom_count);

    auto result =
        WidthDerivatives{.pairs = Eigen::MatrixXd::Zero(n, n), .sources = Eigen::VectorXd::Zero(n)};
    for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
        const auto i = static_cast<Eigen::Index>(atom_index);
        const auto& position = geometry.position(atom_index);
        for (const auto& source : input.fixed_sources()) {
            result.sources(i) +=
                source.charge * width[atom_index] *
                interaction_width_factor(core::squared_distance(position, source.position),
                                         width[atom_index], 0.0);
        }
        for (std::size_t other_atom_index = atom_index + 1; other_atom_index < atom_count;
             ++other_atom_index) {
            const auto j = static_cast<Eigen::Index>(other_atom_index);
            const auto factor =
                interaction_width_factor(geometry.squared_distance(atom_index, other_atom_index),
                                         width[atom_index], width[other_atom_index]);
            result.pairs(i, j) = width[atom_index] * factor;
            result.pairs(j, i) = width[other_atom_index] * factor;
        }
    }
    return result;
}

// The atom channels in the order of the methods' requirements().
constexpr std::size_t electronegativity_channel = 0;
constexpr std::size_t hardness_channel = 1;
constexpr std::size_t width_channel = 2;

class SqeLinearization final : public LinearizedCharges {
  public:
    // An input without bonds has no system.
    SqeLinearization(Eigen::VectorXd charges, std::optional<sqe_core::System> system,
                     WidthDerivatives widths)
        : charges_{std::move(charges)}, system_{std::move(system)}, widths_{std::move(widths)} {
        if (system_.has_value()) {
            // The parts of u that do not depend on the column.
            split_atom_charges_ = charges_;
            if (system_->initial_charges.size() != 0) {
                split_atom_charges_ -= system_->initial_charges;
            }
            width_terms_ = widths_.sources + widths_.pairs * charges_;
        }
    }

    [[nodiscard]] auto charges() const noexcept -> std::span<const double> override {
        return {charges_.data(), static_cast<std::size_t>(charges_.size())};
    }

    [[nodiscard]] auto jacobian(const ParameterColumns& columns) const
        -> std::vector<double> override {
        return jacobian(columns, {});
    }

    // jacobian(), with the initial charges changing too: by column j of initial_charges.values
    // along column initial_charges.columns[j], for variants whose initial charges depend on
    // parameters.
    [[nodiscard]] auto jacobian(const ParameterColumns& columns,
                                const ColumnBlock& initial_charges) const -> std::vector<double> {
        const auto n = charges_.size();
        const auto bond_count = system_.has_value() ? system_->bonds.size() : 0;
        check_columns(columns, 3, 1, 0, static_cast<std::size_t>(n), bond_count, "SQE-family");
        const auto k = static_cast<Eigen::Index>(columns.count);

        // dq starts as d q_initial.
        Eigen::MatrixXd result = Eigen::MatrixXd::Zero(n, k);
        const auto initial_count = initial_charges.columns.size();
        for (std::size_t index = 0; index < initial_count; ++index) {
            result.col(static_cast<Eigen::Index>(initial_charges.columns[index])) +=
                initial_charges.values.col(static_cast<Eigen::Index>(index));
        }

        if (system_.has_value()) {
            const auto& bonds = system_->bonds;
            const auto& p = system_->split_charges;
            const auto m = p.size();

            // u, one column per parameter column.
            Eigen::MatrixXd u = Eigen::MatrixXd::Zero(n, k);
            const auto subtract = [&u](const std::vector<std::size_t>& assignment,
                                       const Eigen::VectorXd& values) {
                for (std::size_t atom = 0; atom < assignment.size(); ++atom) {
                    if (assignment[atom] != ParameterColumns::none) {
                        const auto i = static_cast<Eigen::Index>(atom);
                        u(i, static_cast<Eigen::Index>(assignment[atom])) -= values(i);
                    }
                }
            };
            subtract(columns.atom[electronegativity_channel], Eigen::VectorXd::Ones(n));
            subtract(columns.atom[hardness_channel], split_atom_charges_);
            subtract(columns.atom[width_channel], width_terms_);
            // The dense products only for the columns that have widths or initial charges.
            const auto widths = indicator_block(columns.atom[width_channel]);
            const Eigen::MatrixXd width_products =
                widths_.pairs.transpose() * (charges_.asDiagonal() * widths.values);
            for (std::size_t index = 0; index < widths.columns.size(); ++index) {
                u.col(static_cast<Eigen::Index>(widths.columns[index])) -=
                    width_products.col(static_cast<Eigen::Index>(index));
            }
            if (initial_count != 0) {
                const auto& charge_matrix = system_->charge_matrix;
                const Eigen::MatrixXd initial_products =
                    charge_matrix * initial_charges.values -
                    charge_matrix.diagonal().asDiagonal() * initial_charges.values;
                for (std::size_t index = 0; index < initial_count; ++index) {
                    u.col(static_cast<Eigen::Index>(initial_charges.columns[index])) -=
                        initial_products.col(static_cast<Eigen::Index>(index));
                }
            }

            // S dp = T u - d kappa p.
            Eigen::MatrixXd split = Eigen::MatrixXd(m, k);
            for (Eigen::Index a = 0; a < m; ++a) {
                const auto [i, j] = bonds[static_cast<std::size_t>(a)];
                split.row(a) = u.row(i) - u.row(j);
                const auto column = columns.bond[0][static_cast<std::size_t>(a)];
                if (column != ParameterColumns::none) {
                    split(a, static_cast<Eigen::Index>(column)) -= p(a);
                }
            }
            const auto factor = system_->split_matrix.triangularView<Eigen::Lower>();
            factor.solveInPlace(split);
            factor.transpose().solveInPlace(split);

            // dq += T^T dp.
            for (Eigen::Index a = 0; a < m; ++a) {
                const auto [i, j] = bonds[static_cast<std::size_t>(a)];
                result.row(i) += split.row(a);
                result.row(j) -= split.row(a);
            }
        }

        auto values = std::vector<double>(result.data(), result.data() + result.size());
        require_finite(values, "SQE-family derivative is not finite");
        return values;
    }

  private:
    Eigen::VectorXd charges_;
    std::optional<sqe_core::System> system_;
    WidthDerivatives widths_;
    // T^T p.
    Eigen::VectorXd split_atom_charges_;
    // dphi_k/dwidth_k + sum_j dJ_kj/dwidth_k q_j.
    Eigen::VectorXd width_terms_;
};

// The charges of sqe_core::calculate() with these initial charges, and their derivatives.
[[nodiscard]] auto linearize_sqe(const CalculationInput& input,
                                 const std::span<const double> initial_charge_values)
    -> std::unique_ptr<SqeLinearization> {
    if (input.molecule().bond_count() == 0) {
        // calculate() returns the zero or initial charges.
        const auto charges = sqe_core::calculate(input, initial_charge_values);
        return std::make_unique<SqeLinearization>(
            Eigen::Map<const Eigen::VectorXd>{charges.data(),
                                              static_cast<Eigen::Index>(charges.size())},
            std::nullopt, WidthDerivatives{});
    }

    // The checks and the system of calculate().
    sqe_core::validate_input(input, initial_charge_values);
    auto system = sqe_core::build_system(input, initial_charge_values);
    sqe_core::solve(system);
    Eigen::VectorXd charges = sqe_core::atom_charges(system);
    // The check calculate() gets from AtomicCharges, which rejects non-finite charges.
    static_cast<void>(charges::AtomicCharges{
        std::vector<double>(charges.data(), charges.data() + charges.size())});

    return std::make_unique<SqeLinearization>(std::move(charges), std::move(system),
                                              width_derivatives(input));
}

// The linear part of the SQE+qp normalization of q0 (component_normalized_charges()): subtracts
// each bonded part's mean. The normalization is this projection plus a constant, so a change of q0
// changes the initial charges by the projected change.
auto remove_component_means(const std::vector<std::vector<std::size_t>>& components,
                            const std::span<double> values) -> void {
    for (const auto& component : components) {
        auto sum = 0.0;
        for (const auto atom_index : component) {
            sum += values[atom_index];
        }
        const auto mean = sum / static_cast<double>(component.size());
        for (const auto atom_index : component) {
            values[atom_index] -= mean;
        }
    }
}

// SQE+qp's charges depend on q0 through the initial charges. Explicit reference charges replace
// q0, which then has no effect.
class SqeqpLinearization final : public LinearizedCharges {
  public:
    // `components` are the bonded parts, absent when q0 has no effect.
    SqeqpLinearization(std::unique_ptr<SqeLinearization> sqe,
                       std::optional<std::vector<std::vector<std::size_t>>> components)
        : sqe_{std::move(sqe)}, components_{std::move(components)} {}

    [[nodiscard]] auto charges() const noexcept -> std::span<const double> override {
        return sqe_->charges();
    }

    [[nodiscard]] auto jacobian(const ParameterColumns& columns) const
        -> std::vector<double> override {
        // q0 is the last atom parameter; SQE's columns are the others.
        if (columns.atom.size() != 4) {
            throw std::invalid_argument{
                "SQE+qp jacobian: the columns do not match the method's parameters"};
        }
        auto sqe_columns = columns;
        const auto q0 = std::move(sqe_columns.atom.back());
        sqe_columns.atom.pop_back();
        if (!valid_column_assignment(q0, sqe_->charges().size(), columns.count)) {
            throw std::invalid_argument{
                "SQE+qp jacobian: the columns do not match the method's parameters"};
        }
        if (!components_.has_value()) {
            return sqe_->jacobian(sqe_columns);
        }
        // A change of q0 changes the initial charges by its projection.
        auto initial_charges = indicator_block(q0);
        for (Eigen::Index index = 0; index < initial_charges.values.cols(); ++index) {
            remove_component_means(
                *components_, std::span{initial_charges.values.col(index).data(),
                                        static_cast<std::size_t>(initial_charges.values.rows())});
        }
        return sqe_->jacobian(sqe_columns, initial_charges);
    }

  private:
    std::unique_ptr<SqeLinearization> sqe_;
    std::optional<std::vector<std::vector<std::size_t>>> components_;
};

class SqeDerivatives final : public MethodDerivatives {
  public:
    [[nodiscard]] auto linearize(const CalculationInput& input) const
        -> std::unique_ptr<LinearizedCharges> override {
        return linearize_sqe(input, {});
    }
};

class Sqeq0Derivatives final : public MethodDerivatives {
  public:
    // The formal charges are input, not parameters: their derivative is zero.
    [[nodiscard]] auto linearize(const CalculationInput& input) const
        -> std::unique_ptr<LinearizedCharges> override {
        return linearize_sqe(input, builtin_detail::formal_initial_charges(input.molecule()));
    }
};

class SqeqpDerivatives final : public MethodDerivatives {
  public:
    [[nodiscard]] auto linearize(const CalculationInput& input) const
        -> std::unique_ptr<LinearizedCharges> override {
        const auto initial_charges = sqe_core::parameterized_initial_charges(input);
        auto components = std::optional<std::vector<std::vector<std::size_t>>>{};
        if (initial_charges.from_q0) {
            components = features::connected_components(input.topology());
        }
        return std::make_unique<SqeqpLinearization>(linearize_sqe(input, initial_charges.values),
                                                    std::move(components));
    }
};

} // namespace

auto sqe_derivatives() -> const MethodDerivatives& {
    static const SqeDerivatives derivatives;
    return derivatives;
}

auto sqeq0_derivatives() -> const MethodDerivatives& {
    static const Sqeq0Derivatives derivatives;
    return derivatives;
}

auto sqeqp_derivatives() -> const MethodDerivatives& {
    static const SqeqpDerivatives derivatives;
    return derivatives;
}

} // namespace chargefw::parameterization::detail
