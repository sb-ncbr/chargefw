#include "calculation/cutoff_execution.h"
#include "calculation/fixed_charge_partition.h"
#include "calculation/reduced_execution.h"
#include "calculation/target_execution.h"
#include "methods/applicable_method_execution.h"

#include <chargefw/core/molecule.h>
#include <chargefw/features/conformer_features.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/features/spatial_fragment.h>
#include <chargefw/methods/method.h>
#include <chargefw/parameters/classification/parameter_classification.h>

#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>

namespace chargefw::calculation {
namespace {

[[nodiscard]] auto calculate_target(
    const methods::ApplicableMethod& selected, const features::PreparedMolecule& source,
    const parameters::ParameterClassification* source_classification,
    const std::size_t conformer_index, const double radius, const std::size_t max_threads,
    const detail::ProgressContext& progress_ctx,
    const detail::FixedChargePartitionTarget* partition_target) -> charges::AtomicCharges {
    const auto& source_molecule = source.molecule();
    auto values = std::vector<double>(source_molecule.atom_count());
    const auto charge_context =
        detail::prepare_reduced_charge_context(selected, source, source_classification);
    const auto fixed_sources =
        partition_target == nullptr
            ? std::vector<methods::FixedPointSource>{}
            : detail::materialize_fixed_point_sources(*partition_target, conformer_index);
    const features::ConformerFeatures source_geometry{source_molecule, conformer_index};
    const features::SpatialFragmentBuilder fragment_builder{source, source_geometry};

    ::chargefw::calculation::detail::progress_for_indexed(
        source_molecule.atom_count(), max_threads, progress_ctx,
        [&](const std::size_t center_source_atom_index) {
            const auto failure_message = [&](const std::exception& error) {
                const auto original_atom_index =
                    partition_target == nullptr
                        ? center_source_atom_index
                        : partition_target->active_atom_indices.at(center_source_atom_index);
                return "cutoff fragment around source atom " +
                       std::to_string(original_atom_index + 1) + " failed: " + error.what();
            };
            // Add the fragment to the message and keep the kind of failure; anything else, such
            // as std::bad_alloc, passes unchanged.
            try {
                const auto fragment = fragment_builder.build(center_source_atom_index, radius);
                const auto fragment_charges = detail::calculate_fragment_charges(
                    selected, source_classification, fragment, charge_context, fixed_sources);
                values[center_source_atom_index] =
                    fragment_charges[fragment.center_local_atom_index()];
            } catch (const std::invalid_argument& error) {
                throw std::runtime_error{failure_message(error)};
            } catch (const std::runtime_error& error) {
                throw std::runtime_error{failure_message(error)};
            } catch (const std::logic_error& error) {
                throw std::logic_error{failure_message(error)};
            }
        });

    detail::enforce_conserved_charges(values, charge_context);
    return charges::AtomicCharges{std::move(values)};
}

} // namespace

auto calculate_cutoff_charges(const methods::ApplicableMethod& selected,
                              const features::PreparedMoleculeCollection& molecules,
                              const ExecutionPolicy& policy, const std::size_t max_threads,
                              const CalculationObserver& observer,
                              const detail::FixedChargePartition* fixed_charge_partition)
    -> charges::ChargeSet {
    methods::detail::validate_selected_candidate(selected, molecules,
                                                 fixed_charge_partition != nullptr);
    if (fixed_charge_partition != nullptr) {
        detail::validate_partition_active_molecules(molecules, *fixed_charge_partition);
    }
    methods::detail::validate_coordinate_targets(selected, molecules);
    detail::validate_reduced_request(selected, policy, ExecutionMode::cutoff);
    const auto radius = policy.radius();
    if (!radius.has_value()) {
        throw std::logic_error{"validated cutoff execution policy has no radius"};
    }

    return detail::execute_calculation_targets(
        selected, molecules, ExecutionMode::cutoff, true, detail::ParallelizationLevel::fragments,
        max_threads, observer,
        [&](const features::PreparedMolecule& molecule,
            const parameters::ParameterClassification* classification,
            const std::optional<std::size_t> conformer_index,
            const std::size_t fragment_max_threads, const detail::ProgressContext& target_ctx) {
            const auto* partition_target =
                fixed_charge_partition == nullptr
                    ? nullptr
                    : std::addressof(fixed_charge_partition->targets.at(target_ctx.molecule_index));
            return calculate_target(selected, molecule, classification, *conformer_index, *radius,
                                    fragment_max_threads, target_ctx, partition_target);
        });
}

} // namespace chargefw::calculation
