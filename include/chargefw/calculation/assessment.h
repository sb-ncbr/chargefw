#pragma once

#include <chargefw/calculation/execution_policy.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method_applicability.h>
#include <chargefw/methods/method_requirements.h>
#include <chargefw/parameters/classification/classification_result.h>
#include <chargefw/parameters/models/parameter_set.h>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <string>
#include <unordered_map>
#include <variant>
#include <vector>

namespace chargefw::calculation {

class AssessmentResult;
class CalculationObserver;
class PlanIdentity;
namespace detail {
struct FixedChargePartition;
}
struct AssessmentRequest;
struct ExecutionResult;
// Copies the owned molecule and parameter inputs into the assessment result.
[[nodiscard]] auto assess(const AssessmentRequest& request) -> AssessmentResult;
// Snapshots selection configuration and transfers execution inputs into the assessment result.
// Callers must not inspect the request after passing it to this consuming overload.
[[nodiscard]] auto assess(AssessmentRequest&& request) -> AssessmentResult;

// Immutable concrete execution choice owned by an AssessmentResult. Copies remain tied to the same
// assessment and must not be executed with another assessment's prepared molecules.
class ExecutionPlan {
  public:
    [[nodiscard]] auto method() const noexcept -> const methods::Method&;
    // Null for methods without parameters.
    [[nodiscard]] auto parameter_set() const noexcept -> const parameters::ParameterSet*;
    [[nodiscard]] auto method_options() const noexcept -> const methods::MethodOptions&;
    [[nodiscard]] auto policy() const noexcept -> const ExecutionPolicy&;
    [[nodiscard]] auto warnings() const noexcept -> std::span<const methods::ExecutionIssue>;

  private:
    ExecutionPlan(std::shared_ptr<const PlanIdentity> identity,
                  const methods::ApplicableMethod& candidate, ExecutionPolicy policy,
                  std::vector<methods::ExecutionIssue> warnings);

    friend class AssessmentResult;
    friend auto calculate(const AssessmentResult& assessment, const ExecutionPlan& plan,
                          std::size_t max_threads, const CalculationObserver& observer)
        -> ExecutionResult;
    friend auto estimate(const AssessmentResult& assessment, const ExecutionPlan& plan,
                         std::size_t max_threads) -> methods::ResourceEstimate;

    std::shared_ptr<const PlanIdentity> identity_;
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

// Owns application assessment inputs so adapters and bindings do not need to manage native method
// pointers, parameter spans, or prepared-feature lifetimes.
struct AssessmentRequest {
    core::MoleculeCollection molecules;
    // IDs must be unique across this request, including parameter sets for different methods.
    std::vector<parameters::ParameterSet> parameter_sets;
    std::optional<std::string> method_id;
    // Applies only without method_id; a named method is always assessed.
    MethodScope method_scope = MethodScope::recommended;
    // A parameter-set ID requires method_id so the selected parameterization has an unambiguous
    // scientific model.
    std::optional<std::string> parameter_set_id;
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
    AssessmentResult(AssessmentResult&&) noexcept;
    // Prepared features non-owningly refer to molecules_ or the partition-owned active collection,
    // so replacing these independent owners cannot preserve their lifetime invariant.
    auto operator=(AssessmentResult&&) noexcept -> AssessmentResult& = delete;
    ~AssessmentResult();

    [[nodiscard]] auto applicability_seconds() const noexcept -> double;
    [[nodiscard]] auto molecules() const noexcept -> const core::MoleculeCollection&;

    [[nodiscard]] auto plans() const noexcept -> std::span<const ExecutionPlan>;
    [[nodiscard]] auto rejections() const noexcept -> std::span<const Rejection>;
    [[nodiscard]] auto default_plan() const noexcept -> const ExecutionPlan*;

  private:
    AssessmentResult(core::MoleculeCollection molecules,
                     std::vector<parameters::ParameterSet> parameter_sets,
                     std::unique_ptr<detail::FixedChargePartition> fixed_charge_partition);

    [[nodiscard]] static auto assess_owned(AssessmentRequest request) -> AssessmentResult;

    [[nodiscard]] auto prepared_molecules() const noexcept
        -> const features::PreparedMoleculeCollection&;

    auto
    assess_prepared(std::span<const methods::Method* const> selected_methods,
                    const parameters::ClassificationOptions& classification_options,
                    const ResourcePolicy& resource_policy,
                    const std::unordered_map<std::string, methods::MethodOptions>& method_options,
                    const ExecutionSelection& execution_selection) -> void;

    friend auto assess(const AssessmentRequest& request) -> AssessmentResult;
    friend auto assess(AssessmentRequest&& request) -> AssessmentResult;
    friend auto calculate(const AssessmentResult& assessment, std::size_t max_threads,
                          const CalculationObserver& observer) -> ExecutionResult;
    friend auto calculate(const AssessmentResult& assessment, const ExecutionPlan& plan,
                          std::size_t max_threads, const CalculationObserver& observer)
        -> ExecutionResult;
    friend auto estimate(const AssessmentResult& assessment, const ExecutionPlan& plan,
                         std::size_t max_threads) -> methods::ResourceEstimate;

    std::vector<parameters::ParameterSet> parameter_sets_;
    methods::ApplicabilityResult applicability_;
    std::shared_ptr<const PlanIdentity> plan_identity_;
    std::vector<ExecutionPlan> plans_;
    std::vector<Rejection> rejections_;
    double applicability_seconds_ = 0.0;
    std::unique_ptr<core::MoleculeCollection> molecules_;
    std::unique_ptr<detail::FixedChargePartition> fixed_charge_partition_;
    std::unique_ptr<features::PreparedMoleculeCollection> prepared_molecules_;
};

} // namespace chargefw::calculation
