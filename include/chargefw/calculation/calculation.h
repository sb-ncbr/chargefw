#pragma once

#include <chargefw/calculation/assessment.h>
#include <chargefw/calculation/charge_issues.h>
#include <chargefw/calculation/observer.h>
#include <chargefw/charges/charge_set.h>
#include <chargefw/methods/method_requirements.h>
#include <chargefw/parameters/models/parameter_set.h>

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace chargefw::calculation {

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
    std::vector<ExecutionIssue> warnings;
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
    // Cautions about the calculated charges, in assignment order; empty without charges.
    std::vector<ChargeIssue> charge_warnings;

    [[nodiscard]] auto calculated() const noexcept -> bool {
        return status == ExecutionStatus::success && charges.has_value();
    }
};

// Executes one concrete plan from an assessment without repeating preparation, classification, or
// applicability checks. The plan remains reusable after calculation and may be executed repeatedly.
[[nodiscard]] auto calculate(const ExecutionPlan& plan, std::size_t max_threads = 0,
                             const CalculationObserver& observer = default_calculation_observer())
    -> ExecutionResult;

// Executes `plan` as the overload above, with the values of `candidate` instead of the values of
// the plan's parameter set. The assessment's classifications are reused, so `candidate` must be
// compatible with plan.parameter_set(): the same method, entry keys, and parameter names, in the
// same order (parameters::compatible()). Throws std::invalid_argument if the plan's method has no
// parameters, if `candidate` is not compatible, or if `candidate` has the ID of the plan's
// parameter set. The result reports the candidate's parameter-set ID.
[[nodiscard]] auto calculate(const ExecutionPlan& plan, const parameters::ParameterSet& candidate,
                             std::size_t max_threads = 0,
                             const CalculationObserver& observer = default_calculation_observer())
    -> ExecutionResult;

// Estimates the peak working memory and total operations of executing one plan from an assessment
// with a thread limit, capped at oneTBB's allowed parallelism. Fragments of cutoff and
// cover execution are sized from a typical biomolecular atom density. Memory covers the method
// workspaces that can run concurrently and the retained output charges; it excludes
// already-prepared inputs and is not process RSS.
[[nodiscard]] auto estimate(const ExecutionPlan& plan, std::size_t max_threads = 0)
    -> methods::ResourceEstimate;

// Executes the default plan of an assessment, as calculate(plan) does; without a plan, returns a
// no_executable_plan result.
[[nodiscard]] auto calculate(const AssessmentResult& assessment, std::size_t max_threads = 0,
                             const CalculationObserver& observer = default_calculation_observer())
    -> ExecutionResult;

} // namespace chargefw::calculation
