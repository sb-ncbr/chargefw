#include <chargefw/calculation/calculation.h>

#include "calculation/cover_execution.h"
#include "calculation/cutoff_execution.h"
#include "calculation/fixed_charge_partition.h"
#include "calculation/full_execution.h"
#include "calculation/observer_notifications.h"

#include <chargefw/methods/method.h>

#include <chrono>
#include <cstddef>
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

} // namespace

auto to_string(const ExecutionStatus value) -> std::string_view {
    switch (value) {
    case ExecutionStatus::success:
        return "success";
    case ExecutionStatus::invalid_input_or_request:
        return "invalid_input_or_request";
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

auto calculate(const AssessmentResult& assessment, const ExecutionPlan& plan,
               const std::size_t max_threads, const CalculationObserver& observer)
    -> ExecutionResult {
    if (plan.identity_ != assessment.plan_identity_) {
        throw std::invalid_argument{"execution plan belongs to a different assessment"};
    }

    const auto& selected = *plan.candidate_;
    const auto fixed_charge_partition = assessment.fixed_charge_partition_.get();
    const auto effective = EffectiveCalculation{
        .method_id = std::string{selected.method->id()},
        .parameter_set_id = selected.parameter_set == nullptr
                                ? std::nullopt
                                : std::optional{std::string{selected.parameter_set->id()}},
        .method_options = selected.method_options,
        .execution_policy = plan.policy(),
        .execution_issues = {plan.warnings().begin(), plan.warnings().end()},
        .fixed_ions = fixed_charge_partition == nullptr
                          ? std::nullopt
                          : std::optional{make_fixed_ions_provenance(*fixed_charge_partition)}};

    const auto computation_started = std::chrono::steady_clock::now();
    auto status = ExecutionStatus::success;
    auto calculated_charges = std::optional<charges::ChargeSet>{};
    auto failure_message = std::optional<std::string>{};

    try {
        auto result =
            calculate_request(CalculationRequest{.molecules = assessment.prepared_molecules(),
                                                 .selected = selected,
                                                 .execution_policy = plan.policy(),
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

    const auto computation_seconds =
        std::chrono::duration<double>{std::chrono::steady_clock::now() - computation_started}
            .count();
    return ExecutionResult{
        .status = status,
        .charges = std::move(calculated_charges),
        .rejections = {assessment.rejections_.begin(), assessment.rejections_.end()},
        .effective = effective,
        .failure_message = std::move(failure_message),
        .metrics = {.applicability_seconds = assessment.applicability_seconds_,
                    .computation_seconds = computation_seconds}};
}

auto calculate(const AssessmentResult& assessment, const std::size_t max_threads,
               const CalculationObserver& observer) -> ExecutionResult {
    const auto* plan = assessment.default_plan();
    if (plan == nullptr) {
        return ExecutionResult{
            .status = ExecutionStatus::no_executable_plan,
            .charges = std::nullopt,
            .rejections = {assessment.rejections_.begin(), assessment.rejections_.end()},
            .effective = std::nullopt,
            .failure_message = std::nullopt,
            .metrics = {.applicability_seconds = assessment.applicability_seconds_}};
    }
    return calculate(assessment, *plan, max_threads, observer);
}

} // namespace chargefw::calculation
