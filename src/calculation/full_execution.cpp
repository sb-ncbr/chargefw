#include "calculation/full_execution.h"

#include "calculation/fixed_charge_partition.h"
#include "calculation/target_execution.h"
#include "calculation/target_input.h"
#include "methods/applicable_method_execution.h"

#include <chargefw/core/molecule.h>
#include <chargefw/methods/calculation_input.h>
#include <chargefw/methods/method.h>
#include <chargefw/parameters/classification/parameter_classification.h>

#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

namespace chargefw::calculation {
namespace {

auto validate_calculated_charges(const methods::Method& method,
                                 const features::PreparedMolecule& molecule,
                                 const charges::AtomicCharges& atomic_charges) -> void {
    const auto expected_size = molecule.molecule().atom_count();

    if (atomic_charges.size() != expected_size) {
        throw std::runtime_error{"method '" + std::string{method.id()} + "' produced " +
                                 std::to_string(atomic_charges.size()) + " charges for molecule '" +
                                 std::string{molecule.molecule().name()} + "' with " +
                                 std::to_string(expected_size) + " atoms"};
    }
}

[[nodiscard]] auto calculate_method(const methods::Method& method,
                                    const methods::CalculationInput& input)
    -> charges::AtomicCharges {
    try {
        return method.calculate(input);
    } catch (const std::invalid_argument& error) {
        throw std::runtime_error{error.what()};
    }
}

[[nodiscard]] auto calculate_target(const methods::ApplicableMethod& selected,
                                    const features::PreparedMolecule& molecule,
                                    const parameters::ParameterClassification* classification,
                                    const std::optional<std::size_t> conformer_index,
                                    const detail::FixedChargePartitionTarget* partition_target)
    -> charges::AtomicCharges {
    return detail::with_target_input(
        selected, molecule, classification, conformer_index, partition_target,
        [&](const methods::CalculationInput& input) {
            auto atomic_charges = calculate_method(*selected.method, input);
            validate_calculated_charges(*selected.method, molecule, atomic_charges);
            return atomic_charges;
        });
}

} // namespace

auto calculate_full_charges(const methods::ApplicableMethod& selected,
                            const features::PreparedMoleculeCollection& molecules,
                            const std::size_t max_threads, const CalculationObserver& observer,
                            const detail::FixedChargePartition* fixed_charge_partition)
    -> charges::ChargeSet {
    methods::detail::validate_selected_candidate(selected, molecules,
                                                 fixed_charge_partition != nullptr);
    if (fixed_charge_partition != nullptr) {
        detail::validate_partition_active_molecules(molecules, *fixed_charge_partition);
    }
    methods::detail::validate_coordinate_targets(selected, molecules);

    return detail::execute_calculation_targets(
        selected, molecules, ExecutionMode::full, selected.method->requirements().coordinates,
        detail::ParallelizationLevel::targets, max_threads, observer,
        [&](const features::PreparedMolecule& molecule,
            const parameters::ParameterClassification* classification,
            const std::optional<std::size_t> conformer_index, const std::size_t,
            const detail::ProgressContext& context) {
            const auto* partition_target =
                fixed_charge_partition == nullptr
                    ? nullptr
                    : std::addressof(fixed_charge_partition->targets.at(context.molecule_index));
            return calculate_target(selected, molecule, classification, conformer_index,
                                    partition_target);
        });
}

} // namespace chargefw::calculation
