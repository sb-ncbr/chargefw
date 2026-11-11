#pragma once

#include "calculation/fixed_charge_partition.h"

#include <chargefw/calculation/execution_policy.h>
#include <chargefw/calculation/observer.h>
#include <chargefw/charges/atomic_charges.h>
#include <chargefw/charges/charge_collection.h>
#include <chargefw/core/molecule.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/features/spatial_fragment.h>
#include <chargefw/methods/calculation_input.h>
#include <chargefw/methods/method_applicability.h>
#include <chargefw/methods/method_requirements.h>
#include <chargefw/parameters/classification/parameter_classification.h>

#include <cstddef>
#include <exception>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace chargefw::calculation::detail {

struct ProgressContext;

struct ReducedChargeContext {
    std::vector<double> reference_charges;
    std::vector<std::vector<std::size_t>> conservation_groups;
    // The method starts from the reference charges, so each fragment receives those of its atoms
    // and starts where full execution would.
    bool fragments_start_from_references = false;
};

// Projects whole-molecule parameter entry indices to an induced fragment without reclassification.
[[nodiscard]] auto project_classification(const parameters::ParameterClassification& source,
                                          const features::SpatialFragment& fragment)
    -> parameters::ParameterClassification;
auto validate_reduced_request(const methods::ApplicableMethod& selected,
                              const ExecutionPolicy& policy, ExecutionMode mode) -> void;
[[nodiscard]] auto prepare_reduced_charge_context(
    const methods::ApplicableMethod& selected, const features::PreparedMolecule& source,
    const parameters::ParameterClassification* source_classification) -> ReducedChargeContext;
[[nodiscard]] auto fragment_target_charge(const ReducedChargeContext& context,
                                          const features::SpatialFragment& fragment) -> double;
auto enforce_conserved_charges(std::vector<double>& values, const ReducedChargeContext& context)
    -> void;
// fixed_sources is borrowed and forwarded unchanged during the synchronous method call.
[[nodiscard]] auto calculate_fragment_charges(
    const methods::ApplicableMethod& selected,
    const parameters::ParameterClassification* source_classification,
    const features::SpatialFragment& fragment, const ReducedChargeContext& charge_context,
    std::span<const methods::FixedPointSource> fixed_sources = {}) -> charges::AtomicCharges;

template <typename CalculateFragment>
auto with_fragment_failure_context(const std::string_view label, const std::size_t atom_index,
                                   const FixedChargePartitionTarget* partition_target,
                                   CalculateFragment&& calculate_fragment) -> void {
    const auto failure_message = [&](const std::exception& error) {
        const auto original_atom_index = partition_target == nullptr
                                             ? atom_index
                                             : partition_target->active_atom_indices.at(atom_index);
        return std::string{label} + " fragment around source atom " +
               std::to_string(original_atom_index + 1) + " failed: " + error.what();
    };
    // Method input failures are calculation failures; other logic errors retain their kind.
    // Exceptions outside these categories, such as std::bad_alloc, pass unchanged.
    try {
        calculate_fragment();
    } catch (const std::invalid_argument& error) {
        throw std::runtime_error{failure_message(error)};
    } catch (const std::runtime_error& error) {
        throw std::runtime_error{failure_message(error)};
    } catch (const std::logic_error& error) {
        throw std::logic_error{failure_message(error)};
    }
}

using ReducedTargetCalculator =
    auto (*)(const methods::ApplicableMethod& selected, const features::PreparedMolecule& molecule,
             const parameters::ParameterClassification* classification, std::size_t conformer_index,
             double radius, std::size_t max_threads, const ProgressContext& progress_ctx,
             const FixedChargePartitionTarget* partition_target) -> charges::AtomicCharges;

[[nodiscard]] auto
calculate_reduced_charges(const ExecutionMode mode, const methods::ApplicableMethod& selected,
                          const features::PreparedMoleculeCollection& molecules,
                          const ExecutionPolicy& policy, const std::size_t max_threads,
                          const CalculationObserver& observer,
                          const FixedChargePartition* fixed_charge_partition,
                          ReducedTargetCalculator calculate_target) -> charges::ChargeSet;

} // namespace chargefw::calculation::detail
