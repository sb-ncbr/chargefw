#include "methods/find_applicable_methods.h"

#include <chargefw/methods/method.h>
#include <chargefw/methods/method_applicability.h>
#include <chargefw/methods/parameter_prerequisites.h>

#include <exception>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace chargefw::methods {

namespace {

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

        auto method_result = check_method_prerequisites(*method, request.molecules, method_options);

        if (!method_result) {
            add_rejected(result, method_index, std::nullopt, std::move(method_result.issues));

            continue;
        }

        const auto requirements = method->requirements();

        if (!requirements.requires_parameters()) {
            result.applicable.push_back(ApplicableMethod{.method = method,
                                                         .parameter_set = nullptr,
                                                         .method_options = method_options,
                                                         .classifications = {}});

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

            result.applicable.push_back(
                ApplicableMethod{.method = method,
                                 .parameter_set = &parameter_set,
                                 .method_options = method_options,
                                 .classifications = std::move(parameter_result.classifications)});
        }
    }

    return result;
}

} // namespace chargefw::methods
