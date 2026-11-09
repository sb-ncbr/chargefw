#include "methods/builtin/sqeqp.h"

#include "features/topology_helpers.h"
#include "methods/builtin/sqe.h"
#include "methods/builtin/sqe_derivatives.h"
#include "methods/parameter_columns.h"
#include "methods/reference_charges.h"

#include <chargefw/core/molecule.h>
#include <chargefw/parameters/models/parameter_view.h>

#include <Eigen/Core>

#include <algorithm>
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

// The initial charges: the atom parameter q0, normalized per bonded part.
[[nodiscard]] auto normalized_q0(const CalculationInput& input) -> std::vector<double> {
    const auto& molecule = input.molecule();
    const auto initial_charge_parameter = input.parameters().atom("q0");
    auto initial_charges = std::vector<double>(molecule.atom_count());
    for (std::size_t atom_index = 0; atom_index < molecule.atom_count(); ++atom_index) {
        initial_charges[atom_index] = initial_charge_parameter[atom_index];
    }
    return methods::detail::component_normalized_charges(
        molecule, input.topology(), std::move(initial_charges), input.target_charge());
}

// SQE+qp's charges depend on q0 through the initial charges. The normalization adds a constant to
// q0 and subtracts each bonded part's mean, so q0's gradient is the initial charges' gradient with
// each part's mean subtracted. Explicit reference charges replace q0, which then has no effect.
class SQEqpLinearization final : public LinearizedCharges {
  public:
    // `components` is absent when explicit reference charges replace q0.
    SQEqpLinearization(std::unique_ptr<sqe_core::Linearization> sqe,
                       std::optional<std::vector<std::vector<std::size_t>>> components)
        : sqe_{std::move(sqe)}, components_{std::move(components)} {}

    [[nodiscard]] auto charges() const noexcept -> std::span<const double> override {
        return sqe_->charges();
    }

    [[nodiscard]] auto gradient(const std::span<const double> cotangent) const
        -> ParameterGradient override {
        if (!components_.has_value()) {
            auto result = sqe_->gradient(cotangent);
            result.atom.emplace_back(cotangent.size(), 0.0);
            return result;
        }
        auto [result, q0] = sqe_->gradients(cotangent);
        methods::detail::remove_component_means(*components_, q0);
        if (!std::ranges::all_of(q0, [](const double value) { return std::isfinite(value); })) {
            throw std::runtime_error{"SQE+qp derivative is not finite"};
        }
        result.atom.push_back(std::move(q0));
        return result;
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
        if (q0.size() != sqe_->charges().size() ||
            !std::ranges::all_of(q0, [&columns](const std::size_t column) {
                return column < columns.count || column == ParameterColumns::none;
            })) {
            throw std::invalid_argument{
                "SQE+qp jacobian: the columns do not match the method's parameters"};
        }
        if (!components_.has_value()) {
            return sqe_->jacobian(sqe_columns);
        }
        // A change of q0 changes the initial charges by its projection, as in gradient().
        auto initial_charges = methods::detail::indicator_block(q0);
        for (Eigen::Index index = 0; index < initial_charges.values.cols(); ++index) {
            methods::detail::remove_component_means(
                *components_, std::span{initial_charges.values.col(index).data(),
                                        static_cast<std::size_t>(initial_charges.values.rows())});
        }
        return sqe_->jacobian(sqe_columns, initial_charges);
    }

  private:
    std::unique_ptr<sqe_core::Linearization> sqe_;
    std::optional<std::vector<std::vector<std::size_t>>> components_;
};

} // namespace

auto SQEqpMethod::calculate(const CalculationInput& input) const -> charges::AtomicCharges {
    if (input.molecule().atom_count() == 0) {
        return charges::AtomicCharges{sqe_core::calculate(input)};
    }

    if (const auto references = input.reference_charges(); !references.empty()) {
        return charges::AtomicCharges{sqe_core::calculate(input, references)};
    }

    return charges::AtomicCharges{sqe_core::calculate(input, normalized_q0(input))};
}

auto SQEqpMethod::linearize(const CalculationInput& input) const
    -> std::unique_ptr<LinearizedCharges> {
    if (input.molecule().atom_count() == 0) {
        return std::make_unique<SQEqpLinearization>(sqe_core::linearize(input, {}), std::nullopt);
    }

    if (const auto references = input.reference_charges(); !references.empty()) {
        return std::make_unique<SQEqpLinearization>(sqe_core::linearize(input, references),
                                                    std::nullopt);
    }

    return std::make_unique<SQEqpLinearization>(
        sqe_core::linearize(input, normalized_q0(input)),
        features::connected_components(input.topology().adjacency()));
}

auto SQEqpMethod::estimate(const TargetSize& target, const MethodOptions& /*options*/) const
    -> ResourceEstimate {
    return sqe_core::estimate(target, true);
}

} // namespace chargefw::methods::builtin
