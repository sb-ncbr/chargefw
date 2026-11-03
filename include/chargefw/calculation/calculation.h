#pragma once

#include <chargefw/calculation/assessment.h>
#include <chargefw/calculation/observer.h>
#include <chargefw/charges/charge_collection.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method_requirements.h>

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

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

struct CalculationMetrics {
    double applicability_seconds = 0.0;
    double computation_seconds = 0.0;
};

// The concrete candidate and execution policy selected by the application-facing facade.
struct EffectiveCalculation {
    std::string method_id;
    std::optional<std::string> parameter_set_id;
    methods::MethodOptions method_options;
    ExecutionPolicy execution_policy;
    std::vector<methods::ExecutionIssue> execution_issues;
    std::optional<FixedIons> fixed_ions{};
};

enum class ExecutionStatus : std::uint8_t {
    success,
    no_executable_plan,
    numerical_failure,
    cancelled,
};

[[nodiscard]] auto to_string(ExecutionStatus value) -> std::string_view;

// Result of the application-facing facade, which performs applicability and automatic selection.
struct ExecutionResult {
    ExecutionStatus status = ExecutionStatus::success;
    std::optional<charges::ChargeSet> charges;
    std::vector<Rejection> rejections;
    std::optional<EffectiveCalculation> effective;
    std::optional<std::string> failure_message;
    CalculationMetrics metrics;
    // The input issues of the assessment, in molecule order.
    std::vector<InputIssue> input_warnings;

    [[nodiscard]] auto calculated() const noexcept -> bool {
        return status == ExecutionStatus::success && charges.has_value();
    }

    [[nodiscard]] auto cancelled() const noexcept -> bool {
        return status == ExecutionStatus::cancelled;
    }
};

// Executes an already applicable and selected candidate. Its stored parameter classifications are
// used directly; classification policy is not reconsidered during execution. Emits one enclosing
// computation_started/computation_finished observer pair. Cancellation and calculation failures
// propagate to the caller after computation_finished.
[[nodiscard]] auto calculate(const CalculationRequest& request) -> CalculationResult;

// Executes one concrete plan from an assessment without repeating preparation, classification, or
// applicability checks. The plan must originate from the supplied assessment. Both objects remain
// reusable after calculation and independent plans may be executed repeatedly.
[[nodiscard]] auto calculate(const AssessmentResult& assessment, const ExecutionPlan& plan,
                             std::size_t max_threads = 0,
                             const CalculationObserver& observer = default_calculation_observer())
    -> ExecutionResult;

// Estimates the peak working memory and total operations of executing one plan from an assessment
// with a thread limit, capped at oneTBB's allowed parallelism. Fragments of cutoff and
// cover execution are sized from a typical biomolecular atom density. Memory covers the method
// workspaces that can run concurrently and the retained output charges; it excludes
// already-prepared inputs and is not process RSS.
[[nodiscard]] auto estimate(const AssessmentResult& assessment, const ExecutionPlan& plan,
                            std::size_t max_threads = 0) -> methods::ResourceEstimate;

// Executes a plan returned by assess() without repeating preparation or applicability assessment.
// This convenience overload executes the default plan. It delegates computation observation to the
// reusable plan overload and converts a missing plan to a no_executable_plan result.
[[nodiscard]] auto calculate(const AssessmentResult& assessment, std::size_t max_threads = 0,
                             const CalculationObserver& observer = default_calculation_observer())
    -> ExecutionResult;

} // namespace chargefw::calculation
