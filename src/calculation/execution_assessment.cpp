#include "calculation/execution_assessment.h"

#include <chargefw/methods/method_requirements.h>

#include <cstddef>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace chargefw::calculation {

auto to_string(const ExecutionIssueKind value) -> std::string_view {
    switch (value) {
    case ExecutionIssueKind::resource_threshold_exceeded:
        return "resource_threshold_exceeded";
    case ExecutionIssueKind::unsupported_execution_mode:
        return "unsupported_execution_mode";
    }
    throw std::invalid_argument{"unknown execution issue kind"};
}

namespace detail {

namespace {

[[nodiscard]] auto
has_expensive_full_complexity(const methods::ResourceRequirements& resources) noexcept -> bool {
    using methods::ComplexityTerm;
    const auto cubic = [](const ComplexityTerm complexity) noexcept -> bool {
        return complexity == ComplexityTerm::atoms_cubed ||
               complexity == ComplexityTerm::bonds_cubed ||
               complexity == ComplexityTerm::atoms_plus_bonds_cubed;
    };
    const auto quadratic = [](const ComplexityTerm complexity) noexcept -> bool {
        return complexity == ComplexityTerm::atoms_squared ||
               complexity == ComplexityTerm::bonds_squared ||
               complexity == ComplexityTerm::atoms_plus_bonds_squared ||
               complexity == ComplexityTerm::atoms_squared_plus_bonds_squared;
    };

    return cubic(resources.time) || quadratic(resources.memory);
}

// One issue for every molecule with more atoms than `threshold`; none without a threshold.
[[nodiscard]] auto threshold_issues(const methods::Method& method,
                                    const features::PreparedMoleculeCollection& molecules,
                                    const std::optional<std::size_t> threshold,
                                    const std::string_view execution, const std::string_view limit)
    -> std::vector<ExecutionIssue> {
    auto issues = std::vector<ExecutionIssue>{};
    if (!threshold.has_value()) {
        return issues;
    }
    for (std::size_t molecule_index = 0; molecule_index < molecules.size(); ++molecule_index) {
        if (molecules[molecule_index].molecule().atom_count() > *threshold) {
            issues.push_back(ExecutionIssue{
                .kind = ExecutionIssueKind::resource_threshold_exceeded,
                .message = "method '" + std::string{method.id()} + "' " + std::string{execution} +
                           " execution exceeds the " + std::string{limit} + " threshold of " +
                           std::to_string(*threshold) + " atoms",
                .molecule_index = molecule_index});
        }
    }
    return issues;
}

[[nodiscard]] auto available(std::vector<ExecutionIssue> issues) -> ExecutionAssessment {
    return ExecutionAssessment{.availability = issues.empty()
                                                   ? ExecutionAvailability::available
                                                   : ExecutionAvailability::available_with_warning,
                               .issues = std::move(issues)};
}

} // namespace

auto assess_execution(const methods::Method& method,
                      const features::PreparedMoleculeCollection& molecules,
                      const ResourcePolicy& resource_policy, const ExecutionMode mode)
    -> ExecutionAssessment {
    const auto requirements = method.requirements();
    const auto& resources = requirements.resources;

    if (mode == ExecutionMode::full) {
        return available(has_expensive_full_complexity(resources)
                             ? threshold_issues(method, molecules,
                                                resource_policy.cutoff_atom_threshold, "full",
                                                "cutoff")
                             : std::vector<ExecutionIssue>{});
    }

    const auto cutoff = mode == ExecutionMode::cutoff;
    const auto supported =
        (cutoff ? resources.supports_cutoff : resources.supports_cover) &&
        resources.reduced_charge_policy != methods::ReducedChargePolicy::unsupported;
    if (supported && requirements.coordinates) {
        return available(cutoff ? threshold_issues(method, molecules,
                                                   resource_policy.cover_atom_threshold, "cutoff",
                                                   "cover")
                                : std::vector<ExecutionIssue>{});
    }

    const auto mode_name = std::string{to_string(mode)};
    const auto reason = supported
                            ? " because it does not require coordinates for the radius-based " +
                                  mode_name + " approximation"
                            : "";
    return ExecutionAssessment{
        .availability = ExecutionAvailability::unsupported,
        .issues = {{.kind = ExecutionIssueKind::unsupported_execution_mode,
                    .message = "method '" + std::string{method.id()} + "' does not support " +
                               mode_name + " execution" + reason}}};
}

} // namespace detail

} // namespace chargefw::calculation
