#include <chargefw/methods/method.h>
#include <chargefw/methods/method_applicability.h>
#include <chargefw/methods/parameter_prerequisites.h>

#include <exception>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace chargefw::methods {

auto to_string(const ExecutionAvailability value) -> std::string_view {
    switch (value) {
    case ExecutionAvailability::available:
        return "available";
    case ExecutionAvailability::available_with_warning:
        return "available_with_warning";
    case ExecutionAvailability::unsupported:
        return "unsupported";
    }
    throw std::invalid_argument{"unknown execution availability"};
}

auto to_string(const ExecutionIssueKind value) -> std::string_view {
    switch (value) {
    case ExecutionIssueKind::resource_threshold_exceeded:
        return "resource_threshold_exceeded";
    case ExecutionIssueKind::unsupported_execution_mode:
        return "unsupported_execution_mode";
    }
    throw std::invalid_argument{"unknown execution issue kind"};
}

namespace {

using calculation::ExecutionMode;
using calculation::ResourcePolicy;

[[nodiscard]] auto make_issue(const PrerequisiteIssueKind kind, std::string message)
    -> PrerequisiteIssue {
    return PrerequisiteIssue{.kind = kind, .message = std::move(message)};
}

auto add_rejected(ApplicabilityResult& result, const std::size_t method_index,
                  const std::optional<std::size_t> parameter_set_index,
                  std::vector<PrerequisiteIssue> issues) -> void {
    result.rejected.push_back(RejectedCandidate{.method_index = method_index,
                                                .parameter_set_index = parameter_set_index,
                                                .issues = std::move(issues)});
}

[[nodiscard]] auto copy_issues(std::span<const PrerequisiteIssue> issues)
    -> std::vector<PrerequisiteIssue> {
    return {issues.begin(), issues.end()};
}

[[nodiscard]] auto
valid_classification_count(const MethodRequirements& requirements,
                           const features::PreparedMoleculeCollection& molecules,
                           const CollectionParameterPrerequisiteResult& result) noexcept -> bool {
    if (!requirements.requires_parameter_classification()) {
        return result.classifications.empty();
    }

    return result.classifications.size() == molecules.size();
}

[[nodiscard]] auto make_invalid_classification_count_issue(
    const Method& method, const features::PreparedMoleculeCollection& molecules,
    const CollectionParameterPrerequisiteResult& result) -> PrerequisiteIssue {
    return make_issue(PrerequisiteIssueKind::parameter_classification_failed,
                      "method '" + std::string{method.id()} + "' produced " +
                          std::to_string(result.classifications.size()) +
                          " parameter classifications for " + std::to_string(molecules.size()) +
                          " molecules");
}

[[nodiscard]] auto has_expensive_full_complexity(const ResourceRequirements& resources) noexcept
    -> bool {
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

[[nodiscard]] auto assess_execution(const Method& method,
                                    const features::PreparedMoleculeCollection& molecules,
                                    const ResourcePolicy& resource_policy)
    -> std::vector<ExecutionAssessment> {
    const auto requirements = method.requirements();
    auto full_issues = std::vector<ExecutionIssue>{};

    if (resource_policy.cutoff_atom_threshold.has_value() &&
        has_expensive_full_complexity(requirements.resources)) {
        for (std::size_t molecule_index = 0; molecule_index < molecules.size(); ++molecule_index) {
            const auto atom_count = molecules[molecule_index].molecule().atom_count();
            if (atom_count > *resource_policy.cutoff_atom_threshold) {
                full_issues.push_back(ExecutionIssue{
                    .kind = ExecutionIssueKind::resource_threshold_exceeded,
                    .message = "method '" + std::string{method.id()} +
                               "' full execution exceeds the cutoff threshold of " +
                               std::to_string(*resource_policy.cutoff_atom_threshold) + " atoms",
                    .molecule_index = molecule_index});
            }
        }
    }

    const auto make_reduced_assessment =
        [&method, &requirements, &molecules,
         &resource_policy](const ExecutionMode mode, const bool supported) -> ExecutionAssessment {
        if (supported && requirements.coordinates) {
            auto issues = std::vector<ExecutionIssue>{};
            if (mode == ExecutionMode::cutoff && resource_policy.cover_atom_threshold.has_value()) {
                for (std::size_t molecule_index = 0; molecule_index < molecules.size();
                     ++molecule_index) {
                    if (molecules[molecule_index].molecule().atom_count() >
                        *resource_policy.cover_atom_threshold) {
                        issues.push_back(ExecutionIssue{
                            .kind = ExecutionIssueKind::resource_threshold_exceeded,
                            .message = "method '" + std::string{method.id()} +
                                       "' cutoff execution exceeds the cover threshold of " +
                                       std::to_string(*resource_policy.cover_atom_threshold) +
                                       " atoms",
                            .molecule_index = molecule_index});
                    }
                }
            }
            return ExecutionAssessment{
                .mode = mode,
                .availability = issues.empty() ? ExecutionAvailability::available
                                               : ExecutionAvailability::available_with_warning,
                .issues = std::move(issues)};
        }

        const auto* const mode_name = mode == ExecutionMode::cutoff ? "cutoff" : "cover";
        const auto reason = supported
                                ? " because it does not require coordinates for the radius-based " +
                                      std::string{mode_name} + " approximation"
                                : "";
        return ExecutionAssessment{
            .mode = mode,
            .availability = ExecutionAvailability::unsupported,
            .issues = {{.kind = ExecutionIssueKind::unsupported_execution_mode,
                        .message = "method '" + std::string{method.id()} + "' does not support " +
                                   mode_name + " execution" + reason}}};
    };

    return {ExecutionAssessment{.mode = ExecutionMode::full,
                                .availability = full_issues.empty()
                                                    ? ExecutionAvailability::available
                                                    : ExecutionAvailability::available_with_warning,
                                .issues = std::move(full_issues)},
            make_reduced_assessment(ExecutionMode::cutoff,
                                    requirements.resources.supports_cutoff &&
                                        requirements.resources.reduced_charge_policy !=
                                            ReducedChargePolicy::unsupported),
            make_reduced_assessment(ExecutionMode::cover,
                                    requirements.resources.supports_cover &&
                                        requirements.resources.reduced_charge_policy !=
                                            ReducedChargePolicy::unsupported)};
}

} // namespace

auto find_applicable_methods(const ApplicabilityRequest& request) -> ApplicabilityResult {
    ApplicabilityResult result;

    for (std::size_t method_index = 0; method_index < request.methods.size(); ++method_index) {
        const auto* method = request.methods[method_index];

        if (method == nullptr) {
            throw std::invalid_argument{"applicability check received a null method"};
        }
        // Common parameters belong in a parameter set because they are fitted together with atom
        // or bond parameters; constants on their own are settings.
        if (const auto requirements = method->requirements();
            requirements.requires_common_parameters() &&
            !requirements.requires_parameter_classification()) {
            throw std::invalid_argument{"method '" + std::string{method->id()} +
                                        "' declares only common parameters; declare them as "
                                        "method options instead"};
        }

        MethodOptions method_options;

        try {
            method_options = make_default_options(method->option_schema());
            if (const auto options = request.method_options.find(std::string{method->id()});
                options != request.method_options.end()) {
                for (const auto& [id, value] : options->second.values()) {
                    method_options.set(id, value);
                }
                validate_method_options(method->option_schema(), method_options);
            }
        } catch (const std::exception& error) {
            add_rejected(result, method_index, std::nullopt,
                         {make_issue(PrerequisiteIssueKind::invalid_options, error.what())});

            continue;
        }

        const auto method_result =
            check_method_prerequisites(*method, request.molecules, method_options);

        if (!method_result) {
            add_rejected(result, method_index, std::nullopt, copy_issues(method_result.issues()));

            continue;
        }

        const auto requirements = method->requirements();

        if (!requirements.requires_parameters()) {
            result.applicable.push_back(ApplicableMethod{
                .method = method,
                .parameter_set = nullptr,
                .method_options = method_options,
                .classifications = {},
                .execution_assessments =
                    assess_execution(*method, request.molecules, request.resource_policy)});

            continue;
        }

        if (request.parameter_sets.empty()) {
            add_rejected(
                result, method_index, std::nullopt,
                {make_issue(PrerequisiteIssueKind::missing_parameters,
                            "method '" + std::string{method->id()} +
                                "' requires parameters, but no parameter sets were provided")});

            continue;
        }

        for (std::size_t parameter_set_index = 0;
             parameter_set_index < request.parameter_sets.size(); ++parameter_set_index) {
            const auto& parameter_set = request.parameter_sets[parameter_set_index];

            if (!parameter_set.method_id().empty() && parameter_set.method_id() != method->id()) {
                continue;
            }

            auto parameter_result = check_parameter_prerequisites(
                *method, request.molecules, parameter_set, request.classification_options);

            if (!parameter_result) {
                add_rejected(result, method_index, parameter_set_index,
                             std::move(parameter_result.issues));

                continue;
            }

            if (!valid_classification_count(requirements, request.molecules, parameter_result)) {
                add_rejected(result, method_index, parameter_set_index,
                             {make_invalid_classification_count_issue(*method, request.molecules,
                                                                      parameter_result)});

                continue;
            }

            result.applicable.push_back(ApplicableMethod{
                .method = method,
                .parameter_set = &parameter_set,
                .method_options = method_options,
                .classifications = std::move(parameter_result.classifications),
                .execution_assessments =
                    assess_execution(*method, request.molecules, request.resource_policy)});
        }
    }

    return result;
}

} // namespace chargefw::methods
