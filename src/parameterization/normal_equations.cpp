#include <chargefw/parameterization/normal_equations.h>
#include <chargefw/parameterization/parameter_vector.h>

#include "parameterization/method_derivatives.h"
#include "parameterization/parameter_columns.h"

#include "calculation/assessment_access.h"
#include "calculation/fixed_charge_partition.h"
#include "calculation/parallel_for.h"
#include "calculation/target_execution.h"
#include "calculation/target_input.h"
#include "methods/applicable_method_execution.h"
#include "methods/method_applicability.h"

#include <chargefw/calculation/execution_policy.h>
#include <chargefw/core/molecule.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/methods/calculation_input.h>
#include <chargefw/methods/method.h>
#include <chargefw/methods/method_requirements.h>
#include <chargefw/parameters/classification/parameter_classification.h>

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace chargefw::parameterization {
namespace {

constexpr auto not_finite = "the normal equations are not finite";

using calculation::detail::AssessmentAccess;

constexpr auto no_slot = std::numeric_limits<std::size_t>::max();

// Targets are calculated in parallel in batches of this size; each batch is then added to the
// result in target order, so the arithmetic does not depend on the threads, and only one batch's
// partial sums are kept at a time.
constexpr std::size_t batch_size = 256;

// For each parameter of the method, the layout slot of each entry, or no_slot when the set has no
// such value. Parameters are matched to slots by name.
struct SlotTable {
    std::vector<std::vector<std::size_t>> atom; // [method atom parameter][atom entry]
    std::vector<std::vector<std::size_t>> bond; // [method bond parameter][bond entry]
    std::vector<std::size_t> common;            // [method common parameter]
};

[[nodiscard]] auto make_slot_table(const methods::MethodRequirements& requirements,
                                   const parameters::ParameterSet& set,
                                   const std::span<const ParameterSlot> layout) -> SlotTable {
    const auto channel = [](const std::vector<std::string_view>& names,
                            const std::string_view name) -> std::optional<std::size_t> {
        const auto found = std::ranges::find(names, name);
        if (found == names.end()) {
            return std::nullopt;
        }
        return static_cast<std::size_t>(found - names.begin());
    };

    auto table = SlotTable{
        .atom = std::vector(requirements.atom_parameters.size(),
                            std::vector<std::size_t>(set.atom().size(), no_slot)),
        .bond = std::vector(requirements.bond_parameters.size(),
                            std::vector<std::size_t>(set.bond().size(), no_slot)),
        .common = std::vector<std::size_t>(requirements.common_parameters.size(), no_slot)};
    for (std::size_t slot = 0; slot < layout.size(); ++slot) {
        const auto& [scope, entry, name] = layout[slot];
        switch (scope) {
        case ParameterScope::common:
            if (const auto index = channel(requirements.common_parameters, name)) {
                table.common[*index] = slot;
            }
            break;
        case ParameterScope::atom:
            if (const auto index = channel(requirements.atom_parameters, name)) {
                table.atom[*index][entry] = slot;
            }
            break;
        case ParameterScope::bond:
            if (const auto index = channel(requirements.bond_parameters, name)) {
                table.bond[*index][entry] = slot;
            }
            break;
        }
    }
    return table;
}

// One target's contribution, over only the slots it contains.
struct TargetTerms {
    // The target's rows: its atoms in source order, fixed ions included.
    std::vector<double> charges;
    // The global slots the target contains, ascending.
    std::vector<std::size_t> slots;
    // slots.size() x slots.size(), row-major.
    std::vector<double> gauss_newton;
    std::vector<double> gradient;
    double sum_of_squares = 0.0;
};

[[nodiscard]] auto linearize_method(const detail::MethodDerivatives& derivatives,
                                    const methods::CalculationInput& input)
    -> std::unique_ptr<detail::LinearizedCharges> {
    try {
        return derivatives.linearize(input);
    } catch (const std::invalid_argument& error) {
        // As in ordinary calculation: invalid values for this target are a numerical failure.
        throw std::runtime_error{error.what()};
    }
}

[[nodiscard]] auto target_terms(
    const methods::ApplicableMethod& selected, const detail::MethodDerivatives& derivatives,
    const features::PreparedMolecule& molecule,
    const parameters::ParameterClassification& classification,
    const std::optional<std::size_t> conformer_index,
    const calculation::detail::FixedChargePartitionTarget* partition_target, const SlotTable& table,
    const std::size_t slot_count, const std::span<const double> references,
    const std::span<const double> weights) -> TargetTerms {
    const auto linearized = calculation::detail::with_target_input(
        selected, molecule, &classification, conformer_index, partition_target,
        [&](const methods::CalculationInput& input) {
            return linearize_method(derivatives, input);
        });
    const auto active = linearized->charges();
    const auto atom_count = molecule.molecule().atom_count();
    const auto bond_count = molecule.molecule().bond_count();
    if (active.size() != atom_count) {
        throw std::runtime_error{"method '" + std::string{selected.method->id()} + "' produced " +
                                 std::to_string(active.size()) + " charges for " +
                                 std::to_string(atom_count) + " atoms"};
    }

    // Active atoms are rows at their source positions; fixed ions have their own rows.
    const auto row_of = [partition_target](const std::size_t atom) {
        return partition_target == nullptr ? atom : partition_target->active_atom_indices[atom];
    };
    auto result = TargetTerms{};
    result.charges.assign(references.size(), 0.0);
    for (std::size_t atom = 0; atom < atom_count; ++atom) {
        result.charges[row_of(atom)] = active[atom];
    }
    const auto sources =
        partition_target == nullptr
            ? std::span<const calculation::FixedAtomCharge>{}
            : std::span<const calculation::FixedAtomCharge>{partition_target->sources};
    for (const auto& source : sources) {
        result.charges[source.atom_index] = source.charge;
    }

    // The slots the target contains, and their local positions.
    auto local = std::vector<std::size_t>(slot_count, no_slot);
    const auto add = [&](const std::size_t slot) {
        if (slot != no_slot && local[slot] == no_slot) {
            local[slot] = 0;
            result.slots.push_back(slot);
        }
    };
    for (const auto& entries : table.atom) {
        for (std::size_t atom = 0; atom < atom_count; ++atom) {
            add(entries[classification.atom()[atom]]);
        }
    }
    for (const auto& entries : table.bond) {
        for (std::size_t bond = 0; bond < bond_count; ++bond) {
            add(entries[classification.bond()[bond]]);
        }
    }
    for (const auto slot : table.common) {
        add(slot);
    }
    std::ranges::sort(result.slots);
    for (std::size_t index = 0; index < result.slots.size(); ++index) {
        local[result.slots[index]] = index;
    }
    const auto k = result.slots.size();
    result.gauss_newton.assign(k * k, 0.0);
    result.gradient.assign(k, 0.0);

    // jacobian() with the target's slots as columns gives J, one row per atom; the weighted rows
    // B = W J and residuals t = W r are formed first, so that w^2 is never computed on its own;
    // then A = B^T B and g = B^T t.
    const auto locals = [&](const std::vector<std::size_t>& entries, const auto& classified,
                            const std::size_t count) {
        auto columns = std::vector<std::size_t>(count, detail::ParameterColumns::none);
        for (std::size_t index = 0; index < count; ++index) {
            const auto slot = entries[classified[index]];
            if (slot != no_slot) {
                columns[index] = local[slot];
            }
        }
        return columns;
    };
    auto columns = detail::ParameterColumns{.count = k, .atom = {}, .bond = {}, .common = {}};
    for (const auto& entries : table.atom) {
        columns.atom.push_back(locals(entries, classification.atom(), atom_count));
    }
    for (const auto& entries : table.bond) {
        columns.bond.push_back(locals(entries, classification.bond(), bond_count));
    }
    for (const auto slot : table.common) {
        columns.common.push_back(slot == no_slot ? detail::ParameterColumns::none : local[slot]);
    }
    const auto jacobian = linearized->jacobian(columns);
    if (jacobian.size() != atom_count * k) {
        throw std::logic_error{"method jacobian does not have the shape of its parameters"};
    }
    const Eigen::Map<const Eigen::MatrixXd> charge_derivatives{
        jacobian.data(), static_cast<Eigen::Index>(atom_count), static_cast<Eigen::Index>(k)};

    auto weighted = std::vector<std::size_t>{};
    for (std::size_t atom = 0; atom < atom_count; ++atom) {
        if (weights[row_of(atom)] != 0.0) {
            weighted.push_back(atom);
        }
    }
    const auto rows_count = static_cast<Eigen::Index>(weighted.size());
    Eigen::MatrixXd rows(rows_count, static_cast<Eigen::Index>(k));
    Eigen::VectorXd residuals(rows_count);
    for (Eigen::Index r = 0; r < rows_count; ++r) {
        const auto atom = weighted[static_cast<std::size_t>(r)];
        const auto weight = weights[row_of(atom)];
        residuals(r) = weight * (active[atom] - references[row_of(atom)]);
        result.sum_of_squares += residuals(r) * residuals(r);
        rows.row(r) = weight * charge_derivatives.row(static_cast<Eigen::Index>(atom));
    }
    Eigen::Map<Eigen::MatrixXd>{result.gauss_newton.data(), static_cast<Eigen::Index>(k),
                                static_cast<Eigen::Index>(k)}
        .noalias() = rows.transpose() * rows;
    Eigen::Map<Eigen::VectorXd>{result.gradient.data(), static_cast<Eigen::Index>(k)}.noalias() =
        rows.transpose() * residuals;

    for (const auto& source : sources) {
        const auto weight = weights[source.atom_index];
        if (weight != 0.0) {
            const auto residual = weight * (source.charge - references[source.atom_index]);
            result.sum_of_squares += residual * residual;
        }
    }

    detail::require_finite(result.gauss_newton, not_finite);
    detail::require_finite(result.gradient, not_finite);
    detail::require_finite(std::span{&result.sum_of_squares, 1}, not_finite);
    return result;
}

auto check_rows(const std::span<const double> values, const std::size_t rows,
                const std::string_view name) -> void {
    if (values.size() != rows) {
        throw std::invalid_argument{std::string{name} + " must have " + std::to_string(rows) +
                                    " values, one per charge, not " +
                                    std::to_string(values.size())};
    }
    if (!std::ranges::all_of(values, [](const double value) { return std::isfinite(value); })) {
        throw std::invalid_argument{std::string{name} + " must be finite"};
    }
}

} // namespace

auto has_exact_derivatives(const methods::Method& method) -> bool {
    return detail::derivatives_of(method) != nullptr;
}

auto normal_equations(const calculation::ExecutionPlan& plan,
                      const parameters::ParameterSet& candidate,
                      const std::span<const double> references,
                      const std::span<const double> weights, const std::size_t max_threads)
    -> NormalEquations {
    const auto selected = AssessmentAccess::with_values(plan, candidate);
    if (plan.policy().mode() != calculation::ExecutionMode::full) {
        throw std::invalid_argument{"normal equations need full execution, but the plan uses " +
                                    std::string{to_string(plan.policy().mode())}};
    }
    const auto requirements = selected.method->requirements();
    const auto* derivatives = detail::derivatives_of(*selected.method);
    if (derivatives == nullptr) {
        throw std::invalid_argument{"method '" + std::string{selected.method->id()} +
                                    "' has no exact derivatives"};
    }

    const auto& molecules = AssessmentAccess::molecules(plan);
    const auto* partition = AssessmentAccess::partition(plan);
    methods::detail::validate_selected_candidate(selected, molecules, partition != nullptr);
    if (partition != nullptr) {
        calculation::detail::validate_partition_active_molecules(molecules, *partition);
    }
    methods::detail::validate_coordinate_targets(selected, molecules);

    // Each target's rows follow the previous target's.
    const auto targets =
        calculation::detail::make_calculation_targets(molecules, requirements.coordinates);
    auto offsets = std::vector<std::size_t>{0};
    for (const auto& target : targets) {
        const auto rows =
            partition == nullptr
                ? molecules[target.molecule_index].molecule().atom_count()
                : partition->targets.at(target.molecule_index).active_atom_indices.size() +
                      partition->targets.at(target.molecule_index).sources.size();
        offsets.push_back(offsets.back() + rows);
    }
    check_rows(references, offsets.back(), "references");
    check_rows(weights, offsets.back(), "weights");
    if (std::ranges::any_of(weights, [](const double weight) { return weight < 0.0; })) {
        throw std::invalid_argument{"weights must not be negative"};
    }

    const auto layout = parameter_layout(candidate);
    const auto table = make_slot_table(requirements, candidate, layout);
    const auto slot_count = layout.size();
    auto result = NormalEquations{.charges = std::vector<double>(offsets.back(), 0.0),
                                  .sum_of_squares = 0.0,
                                  .gauss_newton = std::vector<double>(slot_count * slot_count, 0.0),
                                  .gradient = std::vector<double>(slot_count, 0.0)};

    for (std::size_t first = 0; first < targets.size(); first += batch_size) {
        const auto count = std::min(batch_size, targets.size() - first);
        auto terms = std::vector<TargetTerms>(count);
        calculation::detail::parallel_for_indexed(count, max_threads, [&](const std::size_t index) {
            const auto target_index = first + index;
            const auto& target = targets[target_index];
            const auto& molecule = molecules[target.molecule_index];
            const auto rows = offsets[target_index + 1] - offsets[target_index];
            try {
                terms[index] = target_terms(
                    selected, *derivatives, molecule,
                    selected.classifications[target.molecule_index], target.conformer_index,
                    partition == nullptr ? nullptr : &partition->targets.at(target.molecule_index),
                    table, slot_count, references.subspan(offsets[target_index], rows),
                    weights.subspan(offsets[target_index], rows));
            } catch (const std::invalid_argument& error) {
                throw std::invalid_argument{calculation::detail::make_target_failure_message(
                    selected, molecule, target, error)};
            } catch (const std::logic_error& error) {
                throw std::logic_error{calculation::detail::make_target_failure_message(
                    selected, molecule, target, error)};
            } catch (const std::runtime_error& error) {
                throw std::runtime_error{calculation::detail::make_target_failure_message(
                    selected, molecule, target, error)};
            }
        });

        // In target order, whatever the threads.
        for (std::size_t index = 0; index < count; ++index) {
            const auto& part = terms[index];
            std::ranges::copy(part.charges, result.charges.begin() + static_cast<std::ptrdiff_t>(
                                                                         offsets[first + index]));
            result.sum_of_squares += part.sum_of_squares;
            const auto k = part.slots.size();
            for (std::size_t row = 0; row < k; ++row) {
                for (std::size_t column = 0; column < k; ++column) {
                    result.gauss_newton[part.slots[row] * slot_count + part.slots[column]] +=
                        part.gauss_newton[row * k + column];
                }
                result.gradient[part.slots[row]] += part.gradient[row];
            }
        }
    }

    detail::require_finite(result.gauss_newton, not_finite);
    detail::require_finite(result.gradient, not_finite);
    detail::require_finite(std::span{&result.sum_of_squares, 1}, not_finite);
    return result;
}

} // namespace chargefw::parameterization
