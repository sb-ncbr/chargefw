#pragma once

#include "calculation/parallel_for.h"
#include "methods/applicable_method_execution.h"

#include <chargefw/charges/charge_collection.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method.h>
#include <chargefw/methods/method_applicability.h>

#include <chrono>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace chargefw::calculation::detail {

struct CalculationTarget {
    std::size_t molecule_index;
    std::optional<std::size_t> conformer_index;
};

enum class ParallelizationLevel : std::uint8_t {
    targets,
    fragments,
};

[[nodiscard]] inline auto make_target_failure_message(const methods::ApplicableMethod& selected,
                                                      const features::PreparedMolecule& molecule,
                                                      const CalculationTarget& target,
                                                      const std::exception& error) -> std::string {
    auto message = "method '" + std::string{selected.method->id()} + "', molecule " +
                   std::to_string(target.molecule_index + 1) + " ('" +
                   std::string{molecule.molecule().name()} + "')";
    if (target.conformer_index.has_value()) {
        message += ", conformer " + std::to_string(*target.conformer_index + 1);
    }
    return message + " failed: " + error.what();
}

[[nodiscard]] inline auto
make_calculation_targets(const features::PreparedMoleculeCollection& molecules,
                         const bool requires_conformer) -> std::vector<CalculationTarget> {
    auto targets = std::vector<CalculationTarget>{};
    for (std::size_t molecule_index = 0; molecule_index < molecules.size(); ++molecule_index) {
        const auto& molecule = molecules[molecule_index].molecule();
        const auto target_count = requires_conformer ? molecule.conformer_count() : 1;
        for (std::size_t target_index = 0; target_index < target_count; ++target_index) {
            targets.push_back(CalculationTarget{
                .molecule_index = molecule_index,
                .conformer_index =
                    requires_conformer ? std::optional<std::size_t>{target_index} : std::nullopt});
        }
    }
    return targets;
}

// Executes independent molecule/conformer targets and materializes their assignments in input
// order. Full execution parallelizes targets. Reduced execution processes targets serially and
// parallelizes its active target's independent fragments without nested scheduling.
template <typename CalculateTarget>
[[nodiscard]] auto
execute_calculation_targets(const methods::ApplicableMethod& selected,
                            const features::PreparedMoleculeCollection& molecules,
                            const ExecutionMode mode, const bool requires_conformer,
                            const ParallelizationLevel parallelization,
                            const std::size_t max_threads, const CalculationObserver& observer,
                            CalculateTarget&& calculate_target) -> charges::ChargeSet {
    const auto targets = make_calculation_targets(molecules, requires_conformer);
    const auto computation_start = std::chrono::steady_clock::now();
    auto calculated = std::vector<std::optional<charges::ChargeAssignment>>(targets.size());
    const auto target_max_threads =
        parallelization == ParallelizationLevel::targets ? max_threads : 1;
    const auto fragment_max_threads =
        parallelization == ParallelizationLevel::fragments ? max_threads : 1;

    parallel_for_indexed(targets.size(), target_max_threads, [&](const std::size_t target_index) {
        const auto& target = targets[target_index];
        const auto& molecule = molecules[target.molecule_index];
        const auto* classification =
            selected.uses_parameters()
                ? std::addressof(selected.classifications[target.molecule_index])
                : nullptr;
        const ProgressContext target_ctx{
            .observer = observer,
            .mode = mode,
            .method_id = selected.method->id(),
            .target_index = target_index,
            .target_count = targets.size(),
            .molecule_index = target.molecule_index,
            .conformer_index = target.conformer_index,
            .computation_start = computation_start,
        };

        check_cancellation(observer);
        emit_target_event(observer, CalculationPhase::target_started, target_ctx);
        check_cancellation(observer);

        try {
            calculated[target_index] = charges::ChargeAssignment{
                .target = charges::ChargeTarget{.molecule_index = target.molecule_index,
                                                .conformer_index = target.conformer_index},
                .charges = calculate_target(molecule, classification, target.conformer_index,
                                            fragment_max_threads, target_ctx)};
        } catch (const CalculationCancelled&) {
            throw;
        } catch (const std::invalid_argument& error) {
            throw std::invalid_argument{
                make_target_failure_message(selected, molecule, target, error)};
        } catch (const std::logic_error& error) {
            throw std::logic_error{make_target_failure_message(selected, molecule, target, error)};
        } catch (const std::runtime_error& error) {
            // Other exceptions, such as std::bad_alloc, pass unchanged.
            throw std::runtime_error{
                make_target_failure_message(selected, molecule, target, error)};
        }

        emit_target_event(observer, CalculationPhase::target_finished, target_ctx);
    });

    auto assignments = std::vector<charges::ChargeAssignment>{};
    assignments.reserve(calculated.size());
    for (auto& assignment : calculated) {
        assignments.push_back(std::move(*assignment));
    }

    return charges::ChargeSet{std::string{selected.method->id()}, std::move(assignments),
                              methods::detail::parameter_set_id_for(selected)};
}

} // namespace chargefw::calculation::detail
