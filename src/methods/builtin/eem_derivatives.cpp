// Exact derivatives of the EEM charges.
//
// calculate() (eem.cpp) solves M x = b with x = (q, chi_equalized) (eem_system.h). For a weighted
// sum c . q, let
//
//   y = the first n entries of M^-1 (c, 0)      (M is symmetric).
//
// The derivatives are:
//
//   d / d A_k   = -y_k
//   d / d B_k   = -y_k q_k
//   d / d kappa = -(sum_k A_k d/dA_k + sum_k B_k d/dB_k) / kappa = y . (A + B q) / kappa
//
// The last follows from an invariance: dividing every equation by kappa shows that the charges
// depend only on A / kappa and B / kappa, also with fixed sources, so scaling A, B, and kappa
// together changes nothing.

#include "methods/builtin/eem.h"
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

namespace chargefw::methods::builtin {
namespace {

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

    [[nodiscard]] auto gradient(const std::span<const double> cotangent) const
        -> ParameterGradient override {
        if (cotangent.size() != static_cast<std::size_t>(charges_.size())) {
            throw std::invalid_argument{"EEM gradient: the cotangent has the wrong size"};
        }
        for (const auto value : cotangent) {
            if (!std::isfinite(value)) {
                throw std::invalid_argument{"EEM gradient: cotangents must be finite"};
            }
        }

        const auto n = charges_.size();
        if (n == 0) {
            return {.atom = {{}, {}}, .bond = {}, .common = {0.0}};
        }

        Eigen::VectorXd rhs = Eigen::VectorXd::Zero(n + 1);
        rhs.head(n) = Eigen::Map<const Eigen::VectorXd>{cotangent.data(), n};
        const Eigen::VectorXd y = lu_->solve(rhs).head(n);

        auto result = ParameterGradient{.atom = {{}, {}}, .bond = {}, .common = {0.0}};
        result.atom[A_channel] = to_vector(-y);
        result.atom[B_channel] = to_vector(-y.cwiseProduct(charges_));
        result.common[kappa_channel] = y.dot(kappa_terms_);

        for (const auto& channel : result.atom) {
            require_finite(channel);
        }
        require_finite(result.common);
        return result;
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

} // namespace

auto EEMMethod::linearize(const CalculationInput& input) const
    -> std::unique_ptr<LinearizedCharges> {
    detail::validate_fixed_source_input(input, "EEM");
    if (input.molecule().atom_count() == 0) {
        return std::make_unique<EemLinearization>(std::nullopt, input);
    }
    return std::make_unique<EemLinearization>(eem_core::build_system(input), input);
}

} // namespace chargefw::methods::builtin
