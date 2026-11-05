#pragma once

#include <chargefw/calculation/execution_policy.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/methods/method_prerequisites.h>
#include <chargefw/parameters/classification/parameter_classification.h>
#include <chargefw/parameters/models/parameter_set.h>

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace chargefw::methods {

class Method;

enum class ExecutionAvailability : std::uint8_t {
    available,
    available_with_warning,
    unsupported,
};

[[nodiscard]] auto to_string(ExecutionAvailability value) -> std::string_view;

enum class ExecutionIssueKind : std::uint8_t {
    resource_threshold_exceeded,
    unsupported_execution_mode,
};

[[nodiscard]] auto to_string(ExecutionIssueKind value) -> std::string_view;

struct ExecutionIssue {
    ExecutionIssueKind kind;
    std::string message;
    std::optional<std::size_t> molecule_index = std::nullopt;
};

struct ExecutionAssessment {
    calculation::ExecutionMode mode;
    ExecutionAvailability availability;
    std::vector<ExecutionIssue> issues;
};

struct ApplicableMethod {
    const Method* method = nullptr;
    const parameters::ParameterSet* parameter_set = nullptr;

    MethodOptions method_options;
    std::vector<parameters::ParameterClassification> classifications;
    std::vector<ExecutionAssessment> execution_assessments;

    [[nodiscard]] auto uses_parameters() const noexcept -> bool {
        return parameter_set != nullptr;
    }
};

struct RejectedCandidate {
    std::size_t method_index = 0;
    std::optional<std::size_t> parameter_set_index;
    std::vector<PrerequisiteIssue> issues;
};

struct ApplicabilityResult {
    std::vector<ApplicableMethod> applicable;
    std::vector<RejectedCandidate> rejected;
};

} // namespace chargefw::methods
