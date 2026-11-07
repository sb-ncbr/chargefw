// Exact derivatives of the SQE-family charges.
//
// calculate() (sqe.cpp) solves, for atoms i, j and bonds a:
//
//   S p = T r,   q = T^T p + q_initial,   S = T J T^T + diag(kappa),
//   J_ii = hardness_i,   J_ij = I(d_ij, width_i, width_j),
//   r_i = -electronegativity_i - phi_i(width_i) - sum_{j != i} J_ij q_initial_j,
//
// with phi the fixed sources' potential and the last sum only with initial charges. For a weighted
// sum c . q, let
//
//   lambda = S^-1 T c,   mu = T^T lambda,   q_split = T^T p.
//
// Its derivatives are:
//
//   d / d electronegativity_k = -mu_k
//   d / d hardness_k          = -mu_k q_split_k
//   d / d width_k             = -mu_k (dphi_k/dwidth_k + sum_{j != k} dJ_kj/dwidth_k q_j)
//                               - q_k sum_{j != k} dJ_kj/dwidth_k mu_j
//   d / d kappa_a             = -lambda_a p_a
//   d / d q_initial_i         = c_i - sum_{j != i} J_ij mu_j

#include "methods/builtin/sqe_derivatives.h"

#include "methods/builtin/sqe.h"
#include "methods/builtin/sqe_system.h"

#include "methods/fixed_source_validation.h"

#include <chargefw/core/molecule.h>
#include <chargefw/core/position.h>
#include <chargefw/parameters/models/parameter_view.h>

#include <Eigen/Core>

#include <cmath>
#include <cstddef>
#include <memory>
#include <numbers>
#include <optional>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace chargefw::methods::builtin {
namespace {

// interaction() (sqe.cpp) depends on width_i through s = 2 width_i^2 + 2 width_j^2:
// d/dwidth_i erf(d / sqrt(s)) / d = width_i * g with g = -(4 / sqrt(pi)) exp(-d^2 / s) / s^(3/2).
// Returns g, which is the same for both widths; it is zero where interaction() returns 1/r.
[[nodiscard]] auto interaction_width_factor(const double distance, const double width_i,
                                            const double width_j) -> double {
    const auto width_sum = 2.0 * width_i * width_i + 2.0 * width_j * width_j;

    if (distance * distance >= 36.0 * width_sum) {
        return 0.0;
    }

    return -4.0 * std::numbers::inv_sqrtpi * std::exp(-distance * distance / width_sum) /
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
    const auto width = input.parameters().atom("width");
    const auto atom_count = input.molecule().atom_count();
    const auto n = static_cast<Eigen::Index>(atom_count);

    auto result =
        WidthDerivatives{.pairs = Eigen::MatrixXd::Zero(n, n), .sources = Eigen::VectorXd::Zero(n)};
    for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
        const auto i = static_cast<Eigen::Index>(atom_index);
        const auto& position = geometry.position(atom_index);
        for (const auto& source : input.fixed_sources()) {
            const auto distance = core::distance(position, source.position);
            result.sources(i) += source.charge * width[atom_index] *
                                 interaction_width_factor(distance, width[atom_index], 0.0);
        }
        for (std::size_t other_atom_index = atom_index + 1; other_atom_index < atom_count;
             ++other_atom_index) {
            const auto j = static_cast<Eigen::Index>(other_atom_index);
            const auto factor =
                interaction_width_factor(geometry.distance(atom_index, other_atom_index),
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

[[nodiscard]] auto to_vector(const Eigen::VectorXd& values) -> std::vector<double> {
    return {values.data(), values.data() + values.size()};
}

class SqeLinearization final : public sqe_core::Linearization {
  public:
    // An input without bonds has no system.
    SqeLinearization(Eigen::VectorXd charges, std::optional<sqe_core::System> system,
                     WidthDerivatives widths)
        : charges_{std::move(charges)}, system_{std::move(system)}, widths_{std::move(widths)} {
        if (system_.has_value()) {
            // The parts of the gradient that do not depend on the cotangent.
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

    [[nodiscard]] auto gradient(const std::span<const double> cotangent) const
        -> ParameterGradient override {
        return compute(cotangent, false).parameters;
    }

    [[nodiscard]] auto gradients(const std::span<const double> cotangent) const
        -> sqe_core::Gradients override {
        return compute(cotangent, true);
    }

  private:
    [[nodiscard]] auto compute(const std::span<const double> cotangent,
                               const bool initial_charges) const -> sqe_core::Gradients {
        if (cotangent.size() != static_cast<std::size_t>(charges_.size())) {
            throw std::invalid_argument{"SQE-family gradient: the cotangent has the wrong size"};
        }
        for (const auto value : cotangent) {
            if (!std::isfinite(value)) {
                throw std::invalid_argument{"SQE-family gradient: cotangents must be finite"};
            }
        }

        const auto n = charges_.size();
        const Eigen::Map<const Eigen::VectorXd> c{cotangent.data(), n};
        const auto& q = charges_;
        const Eigen::VectorXd zero = Eigen::VectorXd::Zero(n);
        auto result = ParameterGradient{.atom = {to_vector(zero), to_vector(zero), to_vector(zero)},
                                        .bond = {{}},
                                        .common = {}};

        Eigen::VectorXd mu = zero;
        if (system_.has_value()) {
            const auto& bonds = system_->bonds;
            const auto& p = system_->split_charges;
            const auto m = p.size();

            // lambda = S^-1 T c, with the Cholesky factor L of S = L L^T.
            Eigen::VectorXd lambda(m);
            for (Eigen::Index a = 0; a < m; ++a) {
                const auto [i, j] = bonds[static_cast<std::size_t>(a)];
                lambda(a) = c(i) - c(j);
            }
            const auto factor = system_->split_matrix.triangularView<Eigen::Lower>();
            factor.solveInPlace(lambda);
            factor.transpose().solveInPlace(lambda);

            // mu = T^T lambda.
            for (Eigen::Index a = 0; a < m; ++a) {
                const auto [i, j] = bonds[static_cast<std::size_t>(a)];
                mu(i) += lambda(a);
                mu(j) -= lambda(a);
            }

            const Eigen::VectorXd width =
                -mu.cwiseProduct(width_terms_) - q.cwiseProduct(widths_.pairs * mu);
            result.atom[electronegativity_channel] = to_vector(-mu);
            result.atom[hardness_channel] = to_vector(-mu.cwiseProduct(split_atom_charges_));
            result.atom[width_channel] = to_vector(width);
            result.bond[0] = to_vector(-lambda.cwiseProduct(p));
        }

        auto initial_charge_gradient = std::vector<double>{};
        if (initial_charges) {
            Eigen::VectorXd values = c;
            if (system_.has_value()) {
                // J_offdiagonal mu, as build_system() applies J to the initial charges.
                const auto& charge_matrix = system_->charge_matrix;
                values -= charge_matrix * mu - charge_matrix.diagonal().cwiseProduct(mu);
            }
            initial_charge_gradient = to_vector(values);
            require_finite(initial_charge_gradient);
        }

        for (const auto& channel : result.atom) {
            require_finite(channel);
        }
        require_finite(result.bond[0]);
        return {.parameters = std::move(result),
                .initial_charges = std::move(initial_charge_gradient)};
    }

    static auto require_finite(const std::span<const double> values) -> void {
        for (const auto value : values) {
            if (!std::isfinite(value)) {
                throw std::runtime_error{"SQE-family derivative is not finite"};
            }
        }
    }

    Eigen::VectorXd charges_;
    std::optional<sqe_core::System> system_;
    WidthDerivatives widths_;
    // T^T p.
    Eigen::VectorXd split_atom_charges_;
    // dphi_k/dwidth_k + sum_j dJ_kj/dwidth_k q_j.
    Eigen::VectorXd width_terms_;
};

} // namespace

auto sqe_core::linearize(const CalculationInput& input,
                         const std::span<const double> initial_charge_values)
    -> std::unique_ptr<Linearization> {
    if (input.molecule().bond_count() == 0) {
        // calculate() returns the zero or initial charges.
        const auto charges = calculate(input, initial_charge_values);
        return std::make_unique<SqeLinearization>(
            Eigen::Map<const Eigen::VectorXd>{charges.data(),
                                              static_cast<Eigen::Index>(charges.size())},
            std::nullopt, WidthDerivatives{});
    }

    // The checks and the system of calculate().
    methods::detail::validate_fixed_source_input(input, "SQE-family");
    if (!initial_charge_values.empty() &&
        initial_charge_values.size() != input.molecule().atom_count()) {
        throw std::logic_error{"SQE-family initial charge count must match atom count"};
    }
    auto system = build_system(input, initial_charge_values);
    solve(system);
    Eigen::VectorXd charges = atom_charges(system);

    return std::make_unique<SqeLinearization>(std::move(charges), std::move(system),
                                              width_derivatives(input));
}

} // namespace chargefw::methods::builtin
