#include <chargefw/methods/method.h>
#include <chargefw/methods/method_prerequisites.h>

#include "methods/parameter_issue_messages.h"

#include <stdexcept>

namespace chargefw::methods {

auto to_string(const PrerequisiteIssueKind value) -> std::string_view {
    switch (value) {
    case PrerequisiteIssueKind::invalid_options:
        return "invalid_options";
    case PrerequisiteIssueKind::missing_feature:
        return "missing_feature";
    case PrerequisiteIssueKind::invalid_geometry:
        return "invalid_geometry";
    case PrerequisiteIssueKind::unsupported_molecule:
        return "unsupported_molecule";
    case PrerequisiteIssueKind::missing_parameters:
        return "missing_parameters";
    case PrerequisiteIssueKind::unknown_parameters:
        return "unknown_parameters";
    case PrerequisiteIssueKind::parameter_classification_failed:
        return "parameter_classification_failed";
    case PrerequisiteIssueKind::unsupported_fixed_ions:
        return "unsupported_fixed_ions";
    }
    throw std::invalid_argument{"unknown prerequisite issue kind"};
}

auto check_method_prerequisites(const Method& method,
                                const features::PreparedMoleculeCollection& molecules,
                                const MethodOptions& method_options) -> PrerequisiteResult {
    auto result = PrerequisiteResult{};

    for (std::size_t molecule_index = 0; molecule_index < molecules.size(); ++molecule_index) {
        const auto molecule_result = method.check_method_prerequisites(
            {.prepared_molecule = molecules[molecule_index], .method_options = method_options});

        for (const auto& issue : molecule_result.issues) {
            result.issues.push_back(detail::with_molecule_prefix(
                molecules[molecule_index].molecule(), molecule_index, issue));
        }
    }

    return result;
}

} // namespace chargefw::methods
