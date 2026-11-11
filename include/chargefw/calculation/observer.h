#pragma once

#include <chargefw/calculation/execution_policy.h>

#include <cstddef>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string_view>

namespace chargefw::calculation {

// Coarse calculation phases reported by calculation execution and its executors. Target-level
// events are always emitted; fragment-level events are throttled to avoid excessive callback
// overhead on large inputs.
enum class CalculationPhase : std::uint8_t {
    // Execution-level computation phases. computation_finished is the terminal event, emitted once
    // after every emitted computation_started regardless of success, cancellation, or failure.
    computation_started,
    computation_finished,
    // Target tier: one molecule+conformer pair.
    target_started,
    target_finished,
    // Fragment tier: throttled aggregate completion progress for center-atom (cutoff) or pivot
    // (cover) fragments.
    fragment_progress,
};

[[nodiscard]] auto to_string(CalculationPhase value) -> std::string_view;

// Structured progress snapshot. Fields are only populated for the phase that produced the event.
// The non-owning method_id view is valid only during the on_progress callback; callers that need to
// retain it must copy into owned storage.
struct CalculationProgress {
    CalculationPhase phase{};

    // Effective execution mode of the running calculation.
    ExecutionMode mode{};

    // Selected method ID; valid from computation_started onward.
    std::string_view method_id{};

    // Target tier: current and total target indices across all molecules/conformers.
    std::size_t target_index{};
    std::size_t target_count{};

    // Fragment tier: aggregate completed and total fragment counts within the current target.
    // completed_fragment_count is in [1, fragment_count], may skip values due to throttling, and is
    // not a fragment identifier or zero-based index.
    std::size_t completed_fragment_count{};
    std::size_t fragment_count{};

    // Molecule and conformer identifying the current target.
    std::size_t molecule_index{};
    std::optional<std::size_t> conformer_index{};

    // Elapsed computation seconds at the time of the event.
    double elapsed_seconds{};
};

// Cooperative, per-request observer of calculation progress. The observer is non-owning: callers
// retain the object and pass a reference through calculation requests. Executor callbacks can run
// from oneTBB worker threads; implementations must be thread-safe.
//
// The observer is purely observational: it must not mutate method options, parameters, execution
// policy, geometry, or selection. Implementations should not throw; callback exceptions are ignored
// so they cannot affect calculation control flow. Every calculate() overload emits
// computation_finished after computation_started, including when a non-cancellation calculation
// exception propagates. Returning true from cancelled() requests early termination; the application
// facade converts this into a cancelled result with no partial charges.
class CalculationObserver {
  public:
    CalculationObserver() = default;
    CalculationObserver(const CalculationObserver&) = delete;
    auto operator=(const CalculationObserver&) -> CalculationObserver& = delete;
    CalculationObserver(CalculationObserver&&) = delete;
    auto operator=(CalculationObserver&&) -> CalculationObserver& = delete;
    virtual ~CalculationObserver() = default;

    virtual void on_progress(const CalculationProgress& /*progress*/) const {}

    // When true, the running calculation is asked to stop at its next check point: before each
    // molecule or conformer target and before each cutoff or cover fragment. A full solve of one
    // target is not interrupted.
    [[nodiscard]] virtual auto cancelled() const noexcept -> bool {
        return false;
    }
};

// Stateless observer used when callers do not need progress or cancellation handling.
[[nodiscard]] inline auto default_calculation_observer() -> const CalculationObserver& {
    static const auto observer = CalculationObserver{};
    return observer;
}

// Thrown when a cooperative cancellation request is observed mid-calculation. The facade catches
// it and converts the result to a clean cancellation rather than a partial charge set.
class CalculationCancelled : public std::runtime_error {
  public:
    CalculationCancelled() : std::runtime_error{"calculation was cancelled"} {}
};

} // namespace chargefw::calculation
