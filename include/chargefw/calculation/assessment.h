#pragma once

#include <chargefw/calculation/execution_policy.h>
#include <chargefw/calculation/input_issues.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method_applicability.h>
#include <chargefw/methods/method_requirements.h>
#include <chargefw/parameters/classification/classification_result.h>
#include <chargefw/parameters/classification/parameter_classification.h>
#include <chargefw/parameters/models/parameter_set.h>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <unordered_map>
#include <variant>
#include <vector>

namespace chargefw::calculation {

class AssessmentResult;
class CalculationObserver;
namespace detail {
struct AssessmentState;
struct AssessmentAccess;
struct FixedChargePartition;
} // namespace detail
struct AssessmentRequest;
struct ExecutionResult;
// Copies the owned molecule and parameter inputs into the assessment result.
[[nodiscard]] auto assess(const AssessmentRequest& request) -> AssessmentResult;
// Snapshots selection configuration and transfers execution inputs into the assessment result.
// Callers must not inspect the request after passing it to this consuming overload.
[[nodiscard]] auto assess(AssessmentRequest&& request) -> AssessmentResult;

// Immutable concrete execution choice of an AssessmentResult. A plan, and every copy of it, keeps
// the assessed molecules and parameters alive, so it can be calculated after its assessment is
// destroyed.
class ExecutionPlan {
  public:
    [[nodiscard]] auto method() const noexcept -> const methods::Method&;
    // Null for methods without parameters.
    [[nodiscard]] auto parameter_set() const noexcept -> const parameters::ParameterSet*;
    [[nodiscard]] auto method_options() const noexcept -> const methods::MethodOptions&;
    [[nodiscard]] auto policy() const noexcept -> const ExecutionPolicy&;
    [[nodiscard]] auto warnings() const noexcept -> std::span<const methods::ExecutionIssue>;

  private:
    ExecutionPlan(std::shared_ptr<const detail::AssessmentState> state,
                  const methods::ApplicableMethod& candidate, ExecutionPolicy policy,
                  std::vector<methods::ExecutionIssue> warnings);

    // The parameter entries matched during assessment, one classification per calculated molecule.
    // Its atom indices cover every atom when the parameter set has atom entries, and its bond
    // indices every bond when it has bond entries; otherwise that list is empty, as the bond list
    // is for EEM. Empty for a method without parameters. With fixed ions, the molecules are the
    // calculated parts, without the ions.
    [[nodiscard]] auto classifications() const noexcept
        -> std::span<const parameters::ParameterClassification>;

    // The method of the plan with the parameter set `candidate`, which must have the method, keys,
    // and parameter names of the plan's set and another ID; nothing without a candidate. The
    // classifications of the plan stay valid for the candidate.
    [[nodiscard]] auto with_values(const parameters::ParameterSet* candidate) const
        -> std::optional<methods::ApplicableMethod>;

    // Executes the plan with its own parameter set, or with the values of `candidate` when given.
    [[nodiscard]] auto execute(const parameters::ParameterSet* candidate, std::size_t max_threads,
                               const CalculationObserver& observer) const -> ExecutionResult;

    friend class AssessmentResult;
    friend struct detail::AssessmentAccess;
    friend auto calculate(const ExecutionPlan& plan, std::size_t max_threads,
                          const CalculationObserver& observer) -> ExecutionResult;
    friend auto calculate(const ExecutionPlan& plan, const parameters::ParameterSet& candidate,
                          std::size_t max_threads, const CalculationObserver& observer)
        -> ExecutionResult;
    friend auto estimate(const ExecutionPlan& plan, std::size_t max_threads)
        -> methods::ResourceEstimate;

    std::shared_ptr<const detail::AssessmentState> state_;
    // Points into state_.
    const methods::ApplicableMethod* candidate_ = nullptr;
    ExecutionPolicy policy_{};
    std::vector<methods::ExecutionIssue> warnings_;
};

using RejectionIssue = std::variant<methods::PrerequisiteIssue, methods::ExecutionIssue>;

// A scientific candidate or concrete execution policy excluded during assessment. A missing policy
// means the candidate failed before a runnable plan could be formed.
struct Rejection {
    std::string method_id;
    std::optional<std::string> parameter_set_id;
    std::optional<ExecutionPolicy> policy;
    std::vector<RejectionIssue> issues;
};

struct FixedAtomCharge {
    std::size_t molecule_index; // Zero-based in the supplied molecule collection.
    std::size_t atom_index;     // Zero-based in the supplied molecule's atom collection.
    double charge;
};

struct FixedIons {
    // Explicit prescribed charges; native callers are not restricted to the adapter ion catalog.
    std::vector<FixedAtomCharge> sources;
};

// The methods automatic selection may choose from: only recommended methods, or all methods.
enum class MethodScope : std::uint8_t {
    recommended,
    all,
};

[[nodiscard]] auto method_scope_from_string(std::string_view value) -> MethodScope;
[[nodiscard]] auto to_string(MethodScope value) -> std::string_view;

// Owns application assessment inputs so adapters and bindings do not need to manage native method
// pointers, parameter spans, or prepared-feature lifetimes.
struct AssessmentRequest {
    core::MoleculeCollection molecules;
    // IDs must be unique across this request, including parameter sets for different methods.
    std::vector<parameters::ParameterSet> parameter_sets;
    std::optional<std::string> method_id;
    // A parameter-set ID requires method_id so the selected parameterization has an unambiguous
    // scientific model.
    std::optional<std::string> parameter_set_id;
    // Applies only without method_id; a named method is always assessed.
    MethodScope method_scope = MethodScope::recommended;
    std::unordered_map<std::string, methods::MethodOptions> method_options;
    parameters::ClassificationOptions classification_options{};
    ExecutionSelection execution_selection{};
    ResourcePolicy resource_policy{};
    std::optional<FixedIons> fixed_ions{};
};

// Owns prepared application inputs and reusable concrete execution plans.
class AssessmentResult {
  public:
    AssessmentResult(const AssessmentResult&) = delete;
    auto operator=(const AssessmentResult&) -> AssessmentResult& = delete;
    // The moved-from assessment keeps its molecules, rejections, and warnings, and has no plans.
    AssessmentResult(AssessmentResult&&) noexcept;
    // Not assignable, so that pointers to its plans stay valid for as long as it exists.
    auto operator=(AssessmentResult&&) noexcept -> AssessmentResult& = delete;
    ~AssessmentResult();

    [[nodiscard]] auto applicability_seconds() const noexcept -> double;
    [[nodiscard]] auto molecules() const noexcept -> const core::MoleculeCollection&;

    // Not available on a temporary assessment, whose plans would be destroyed with it; copy a plan
    // to keep it.
    [[nodiscard]] auto plans() const& noexcept -> std::span<const ExecutionPlan>;
    auto plans() const&& -> std::span<const ExecutionPlan> = delete;
    [[nodiscard]] auto rejections() const noexcept -> std::span<const Rejection>;
    // Every input issue of the molecules, whichever plan is executed.
    [[nodiscard]] auto input_warnings() const noexcept -> std::span<const InputIssue>;
    [[nodiscard]] auto default_plan() const& noexcept -> const ExecutionPlan*;
    auto default_plan() const&& -> const ExecutionPlan* = delete;

  private:
    AssessmentResult(std::shared_ptr<const detail::AssessmentState> state,
                     std::vector<ExecutionPlan> plans);

    [[nodiscard]] static auto assess_owned(AssessmentRequest request) -> AssessmentResult;

    // Finds the applicable methods of `state` and records its rejections; returns the plans.
    [[nodiscard]] static auto
    assess_prepared(const std::shared_ptr<detail::AssessmentState>& state,
                    std::span<const methods::Method* const> selected_methods,
                    const parameters::ClassificationOptions& classification_options,
                    const ResourcePolicy& resource_policy,
                    const std::unordered_map<std::string, methods::MethodOptions>& method_options,
                    const ExecutionSelection& execution_selection) -> std::vector<ExecutionPlan>;

    friend auto assess(const AssessmentRequest& request) -> AssessmentResult;
    friend auto assess(AssessmentRequest&& request) -> AssessmentResult;
    friend auto calculate(const AssessmentResult& assessment, std::size_t max_threads,
                          const CalculationObserver& observer) -> ExecutionResult;

    // Never null: a moved-from assessment keeps the state and has no plans.
    std::shared_ptr<const detail::AssessmentState> state_;
    std::vector<ExecutionPlan> plans_;
};

} // namespace chargefw::calculation
