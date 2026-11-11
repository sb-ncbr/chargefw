#include "calculation/cutoff_execution.h"
#include "calculation/fixed_charge_partition.h"
#include "calculation/parallel_for.h"
#include "calculation/reduced_execution.h"

#include <chargefw/core/molecule.h>
#include <chargefw/features/conformer_features.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/features/spatial_fragment.h>
#include <chargefw/methods/method.h>
#include <chargefw/parameters/classification/parameter_classification.h>

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
            detail::with_fragment_failure_context(
                "cutoff", center_source_atom_index, partition_target, [&] {
                    const auto fragment = fragment_builder.build(center_source_atom_index, radius);
                    const auto fragment_charges = detail::calculate_fragment_charges(
                        selected, source_classification, fragment, charge_context, fixed_sources);
                    values[center_source_atom_index] =
                        fragment_charges[fragment.center_local_atom_index()];
                });
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
    return detail::calculate_reduced_charges(ExecutionMode::cutoff, selected, molecules, policy,
                                             max_threads, observer, fixed_charge_partition,
                                             calculate_target);
}

} // namespace chargefw::calculation
