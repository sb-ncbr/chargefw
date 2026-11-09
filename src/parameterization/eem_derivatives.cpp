// Exact derivatives of the EEM charges.
//
// calculate() (eem.cpp) solves M x = b with x = (q, chi_equalized) (eem_system.h). A change of
// the parameters changes the charges by the first n entries of
//
//   M^-1 (-dA - dB q + dkappa (A + B q) / kappa, 0).
//
// The kappa term follows from an invariance: dividing every equation by kappa shows that the
// charges depend only on A / kappa and B / kappa, also with fixed sources, so scaling A, B, and
// kappa together changes nothing. jacobian() takes one such change per column.

#include "parameterization/method_derivatives.h"
#include "parameterization/parameter_columns.h"

#include "methods/builtin/eem_system.h"
#include "methods/fixed_source_validation.h"

#include <chargefw/charges/atomic_charges.h>
#include <chargefw/core/molecule.h>
#include <chargefw/parameters/models/parameter_view.h>

#include <Eigen/Core>
#include <Eigen/LU>

#include <cmath>
#include <cstddef>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace chargefw::parameterization::detail {
namespace {

using methods::CalculationInput;
namespace eem_core = methods::builtin::eem_core;

[[nodiscard]] auto to_vector(const Eigen::VectorXd& values) -> std::vector<double> {
    return {values.data(), values.data() + values.size()};
}

// The atom channels in the order of requirements().
constexpr std::size_t A_channel = 0;
constexpr std::size_t B_channel = 1;
constexpr std::size_t kappa_channel = 0;

class EemLinearization final : public LinearizedCharges {
  public:
    // An input without atoms has no system.
    EemLinearization(std::optional<eem_core::System> system, const CalculationInput& input)
        : matrix_{system.has_value() ? std::move(system->matrix) : Eigen::MatrixXd{}} {
        if (!system.has_value()) {
            return;
        }
        // In place, exactly as solve_in_place() does in calculate().
        lu_.emplace(matrix_);
        const auto n = matrix_.rows() - 1;
        charges_ = lu_->solve(system->rhs).head(n);
        // The check calculate() gets from AtomicCharges, which rejects non-finite charges, for
        // example from a singular system.
        static_cast<void>(charges::AtomicCharges{to_vector(charges_)});

        const auto& parameters = input.parameters();
        const auto A = parameters.atom("A");
        const auto B = parameters.atom("B");
        kappa_terms_.resize(n);
        for (Eigen::Index i = 0; i < n; ++i) {
            const auto atom = static_cast<std::size_t>(i);
            kappa_terms_(i) = A[atom] + B[atom] * charges_(i);
        }
        kappa_terms_ /= parameters.common("kappa");
    }

    [[nodiscard]] auto charges() const noexcept -> std::span<const double> override {
        return {charges_.data(), static_cast<std::size_t>(charges_.size())};
    }

    [[nodiscard]] auto jacobian(const ParameterColumns& columns) const
        -> std::vector<double> override {
        const auto n = charges_.size();
        check_columns(columns, 2, 0, 1, static_cast<std::size_t>(n), 0, "EEM");
        const auto k = static_cast<Eigen::Index>(columns.count);
        if (n == 0) {
            return {};
        }

        Eigen::MatrixXd rhs = Eigen::MatrixXd::Zero(n + 1, k);
        for (Eigen::Index i = 0; i < n; ++i) {
            const auto atom = static_cast<std::size_t>(i);
            if (const auto column = columns.atom[A_channel][atom];
                column != ParameterColumns::none) {
                rhs(i, static_cast<Eigen::Index>(column)) -= 1.0;
            }
            if (const auto column = columns.atom[B_channel][atom];
                column != ParameterColumns::none) {
                rhs(i, static_cast<Eigen::Index>(column)) -= charges_(i);
            }
        }
        if (const auto column = columns.common[kappa_channel]; column != ParameterColumns::none) {
            rhs.col(static_cast<Eigen::Index>(column)).head(n) += kappa_terms_;
        }
        const Eigen::MatrixXd result = lu_->solve(rhs).topRows(n);

        auto values = std::vector<double>(result.data(), result.data() + result.size());
        require_finite(values);
        return values;
    }

  private:
    static auto require_finite(const std::span<const double> values) -> void {
        for (const auto value : values) {
            if (!std::isfinite(value)) {
                throw std::runtime_error{"EEM derivative is not finite"};
            }
        }
    }

    // The LU factors of M after construction; lu_ refers to it, so this object cannot move.
    Eigen::MatrixXd matrix_;
    std::optional<Eigen::PartialPivLU<Eigen::Ref<Eigen::MatrixXd>>> lu_;
    Eigen::VectorXd charges_;
    // (A + B q) / kappa.
    Eigen::VectorXd kappa_terms_;
};

class EemDerivatives final : public MethodDerivatives {
  public:
    [[nodiscard]] auto linearize(const CalculationInput& input) const
        -> std::unique_ptr<LinearizedCharges> override {
        methods::detail::validate_fixed_source_input(input, "EEM");
        if (input.molecule().atom_count() == 0) {
            return std::make_unique<EemLinearization>(std::nullopt, input);
        }
        return std::make_unique<EemLinearization>(eem_core::build_system(input), input);
    }
};

} // namespace

auto eem_derivatives() -> const MethodDerivatives& {
    static const EemDerivatives derivatives;
    return derivatives;
}

} // namespace chargefw::parameterization::detail
