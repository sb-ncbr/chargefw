#include <chargefw/calculation/calculation.h>

#include "calculation/assessment_state.h"
#include "calculation/calculation_request.h"
#include "calculation/charge_checks.h"
#include "calculation/cover_execution.h"
#include "calculation/cutoff_execution.h"
#include "calculation/fixed_charge_partition.h"
#include "calculation/full_execution.h"
#include "calculation/input_checks.h"
#include "calculation/observer_notifications.h"
#include "methods/applicable_method_execution.h"

#include <chargefw/methods/method.h>
#include <chargefw/methods/method_applicability.h>
#include <chargefw/parameters/models/parameter_set.h>

#include <chrono>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace chargefw::calculation {

namespace {

class ComputationFinishedEmitter {
  public:
    ComputationFinishedEmitter(const CalculationObserver& observer, const ExecutionMode mode,
                               const std::string_view method_id,
                               const std::chrono::steady_clock::time_point started) noexcept
        : observer_{observer}, mode_{mode}, method_id_{method_id}, started_{started} {}

    ComputationFinishedEmitter(const ComputationFinishedEmitter&) = delete;
    auto operator=(const ComputationFinishedEmitter&) -> ComputationFinishedEmitter& = delete;
    ComputationFinishedEmitter(ComputationFinishedEmitter&&) = delete;
    auto operator=(ComputationFinishedEmitter&&) -> ComputationFinishedEmitter& = delete;

    ~ComputationFinishedEmitter() noexcept {
        detail::report_progress(observer_, CalculationProgress{
                                               .phase = CalculationPhase::computation_finished,
                                               .mode = mode_,
                                               .method_id = method_id_,
                                               .elapsed_seconds =
                                                   std::chrono::duration<double>{
                                                       std::chrono::steady_clock::now() - started_}
                                                       .count(),
                                           });
    }

  private:
    const CalculationObserver& observer_;
    ExecutionMode mode_;
    std::string_view method_id_;
    std::chrono::steady_clock::time_point started_;
};

[[nodiscard]] auto make_fixed_ions_provenance(const detail::FixedChargePartition& partition)
    -> FixedIons {
    auto provenance = FixedIons{};
    auto source_count = std::size_t{0};
    for (const auto& target : partition.targets) {
        source_count += target.sources.size();
    }
    provenance.sources.reserve(source_count);
    for (const auto& target : partition.targets) {
        provenance.sources.insert(provenance.sources.end(), target.sources.begin(),
                                  target.sources.end());
    }
    return provenance;
}

[[nodiscard]] auto
calculate_request(const CalculationRequest& request,
                  const detail::FixedChargePartition* fixed_charge_partition = nullptr)
    -> CalculationResult {
    if (request.selected.method == nullptr) {
        throw std::invalid_argument{"calculation request has no selected method"};
    }

    const auto computation_started = std::chrono::steady_clock::now();
    const auto computation_finished =
        ComputationFinishedEmitter{request.observer, request.execution_policy.mode(),
                                   request.selected.method->id(), computation_started};
    detail::report_progress(request.observer, CalculationProgress{
                                                  .phase = CalculationPhase::computation_started,
                                                  .mode = request.execution_policy.mode(),
                                                  .method_id = request.selected.method->id(),
                                              });

    switch (request.execution_policy.mode()) {
    case ExecutionMode::full:
        return CalculationResult{.charges = calculate_full_charges(
                                     request.selected, request.molecules, request.max_threads,
                                     request.observer, fixed_charge_partition)};
    case ExecutionMode::cutoff:
        return CalculationResult{
            .charges = calculate_cutoff_charges(request.selected, request.molecules,
                                                request.execution_policy, request.max_threads,
                                                request.observer, fixed_charge_partition)};
    case ExecutionMode::cover:
        return CalculationResult{
            .charges = calculate_cover_charges(request.selected, request.molecules,
                                               request.execution_policy, request.max_threads,
                                               request.observer, fixed_charge_partition)};
    }

    throw std::invalid_argument{"unknown execution policy"};
}

// Close atoms matter only to a method using coordinates. A cancelled calculation reports only the
// assessment warnings, so that it does not start a search.
auto execution_input_warnings(const detail::AssessmentState& state, const methods::Method& method,
                              const ExecutionStatus status) -> std::vector<InputIssue> {
    if (!method.requirements().coordinates || status == ExecutionStatus::cancelled) {
        return state.input_warnings;
    }
    return detail::add_close_atom_issues(state.input_warnings, state.molecules,
                                         state.prepared_molecules,
                                         state.fixed_charge_partition.get());
}

} // namespace

auto to_string(const ExecutionStatus value) -> std::string_view {
    switch (value) {
    case ExecutionStatus::success:
        return "success";
    case ExecutionStatus::no_executable_plan:
        return "no_executable_plan";
    case ExecutionStatus::numerical_failure:
        return "numerical_failure";
    case ExecutionStatus::cancelled:
        return "cancelled";
    }
    throw std::invalid_argument{"unknown execution status"};
}

auto calculate(const CalculationRequest& request) -> CalculationResult {
    return calculate_request(request);
}

auto ExecutionPlan::with_values(const parameters::ParameterSet* candidate) const
    -> std::optional<methods::ApplicableMethod> {
    if (candidate == nullptr) {
        return std::nullopt;
    }
    const auto* assessed = parameter_set();
    if (assessed == nullptr) {
        throw std::invalid_argument{"method '" + std::string{method().id()} +
                                    "' has no parameters to replace"};
    }
    if (!parameters::compatible(*assessed, *candidate)) {
        throw std::invalid_argument{
            "parameter set '" + std::string{candidate->id()} +
            "' does not have the method, keys, and parameter names of parameter set '" +
            std::string{assessed->id()} + "'"};
    }
    if (candidate->id() == assessed->id()) {
        throw std::invalid_argument{"parameter set with changed values must not keep the ID '" +
                                    std::string{assessed->id()} + "'"};
    }
    // The executors, the charge set, and the provenance all read the candidate through this copy.
    auto replaced = *candidate_;
    replaced.parameter_set = candidate;
    return replaced;
}

auto ExecutionPlan::execute(const parameters::ParameterSet* candidate,
                            const std::size_t max_threads,
                            const CalculationObserver& observer) const -> ExecutionResult {
    const auto with_candidate = with_values(candidate);

    const auto& state = *state_;
    const auto& selected = with_candidate.has_value() ? *with_candidate : *candidate_;
    const auto* fixed_charge_partition = state.fixed_charge_partition.get();
    const auto effective = EffectiveCalculation{
        .method_id = std::string{selected.method->id()},
        .parameter_set_id = methods::detail::parameter_set_id_for(selected),
        .method_options = selected.method_options,
        .execution_policy = policy_,
        .execution_issues = warnings_,
        .fixed_ions = fixed_charge_partition == nullptr
                          ? std::nullopt
                          : std::optional{make_fixed_ions_provenance(*fixed_charge_partition)}};

    const auto computation_started = std::chrono::steady_clock::now();
    auto status = ExecutionStatus::success;
    auto calculated_charges = std::optional<charges::ChargeSet>{};
    auto failure_message = std::optional<std::string>{};

    try {
        auto result = calculate_request(CalculationRequest{.molecules = state.prepared_molecules,
                                                           .selected = selected,
                                                           .execution_policy = policy_,
                                                           .max_threads = max_threads,
                                                           .observer = observer},
                                        fixed_charge_partition);
        if (fixed_charge_partition != nullptr) {
            result.charges =
                detail::reassemble_fixed_charge_results(result.charges, *fixed_charge_partition);
        }
        calculated_charges.emplace(std::move(result.charges));
    } catch (const CalculationCancelled&) {
        status = ExecutionStatus::cancelled;
    } catch (const std::invalid_argument&) {
        throw;
    } catch (const std::logic_error& error) {
        // Methods report invalid parameter values and degenerate equations as logic errors.
        status = ExecutionStatus::numerical_failure;
        failure_message = error.what();
    } catch (const std::runtime_error& error) {
        status = ExecutionStatus::numerical_failure;
        failure_message = error.what();
    }
    // Anything else, such as std::bad_alloc, is not a failure of the calculation and propagates.

    auto charge_warnings = calculated_charges.has_value()
                               ? detail::find_charge_issues(state.molecules, *calculated_charges,
                                                            fixed_charge_partition)
                               : std::vector<ChargeIssue>{};
    // Searching molecules with fixed ions anew is part of the computation time.
    auto input_warnings = execution_input_warnings(state, *selected.method, status);
    const auto computation_seconds =
        std::chrono::duration<double>{std::chrono::steady_clock::now() - computation_started}
            .count();
    return ExecutionResult{.status = status,
                           .charges = std::move(calculated_charges),
                           .rejections = state.rejections,
                           .effective = effective,
                           .failure_message = std::move(failure_message),
                           .metrics = {.applicability_seconds = state.applicability_seconds,
                                       .computation_seconds = computation_seconds},
                           .input_warnings = std::move(input_warnings),
                           .charge_warnings = std::move(charge_warnings)};
}

auto calculate(const ExecutionPlan& plan, const std::size_t max_threads,
               const CalculationObserver& observer) -> ExecutionResult {
    return plan.execute(nullptr, max_threads, observer);
}

auto calculate(const ExecutionPlan& plan, const parameters::ParameterSet& candidate,
               const std::size_t max_threads, const CalculationObserver& observer)
    -> ExecutionResult {
    return plan.execute(&candidate, max_threads, observer);
}

auto calculate(const AssessmentResult& assessment, const std::size_t max_threads,
               const CalculationObserver& observer) -> ExecutionResult {
    const auto* plan = assessment.default_plan();
    if (plan == nullptr) {
        const auto& state = *assessment.state_;
        return ExecutionResult{.status = ExecutionStatus::no_executable_plan,
                               .charges = std::nullopt,
                               .rejections = state.rejections,
                               .effective = std::nullopt,
                               .failure_message = std::nullopt,
                               .metrics = {.applicability_seconds = state.applicability_seconds},
                               .input_warnings = state.input_warnings,
                               .charge_warnings = {}};
    }
    return calculate(*plan, max_threads, observer);
}

} // namespace chargefw::calculation
