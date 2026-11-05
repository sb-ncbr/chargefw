#pragma once

#include <chargefw/calculation/execution_policy.h>
#include <chargefw/calculation/observer.h>
#include <chargefw/charges/charge_collection.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method_applicability.h>

#include <cstddef>

namespace chargefw::calculation {

struct CalculationRequest {
    // Non-owning input; the referenced collection and selected candidate must outlive calculate().
    // NOLINTNEXTLINE(cppcoreguidelines-avoid-const-or-ref-data-members)
    const features::PreparedMoleculeCollection& molecules;
    // NOLINTNEXTLINE(cppcoreguidelines-avoid-const-or-ref-data-members)
    const methods::ApplicableMethod& selected;
    ExecutionPolicy execution_policy{};
    // Zero delegates the worker count to the oneTBB runtime.
    std::size_t max_threads = 0;
    // Non-owning observer for progress and cancellation.
    // NOLINTNEXTLINE(cppcoreguidelines-avoid-const-or-ref-data-members)
    const CalculationObserver& observer = default_calculation_observer();
};

struct CalculationResult {
    charges::ChargeSet charges;
};

// Executes an already applicable and selected candidate. Its stored parameter classifications are
// used directly; classification policy is not reconsidered during execution. Emits one enclosing
// computation_started/computation_finished observer pair. Cancellation and calculation failures
// propagate to the caller after computation_finished.
[[nodiscard]] auto calculate(const CalculationRequest& request) -> CalculationResult;

} // namespace chargefw::calculation
