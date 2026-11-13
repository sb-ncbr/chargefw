#pragma once

#include <chargefw/methods/method_options.h>
#include <chargefw/methods/method_prerequisites.h>
#include <chargefw/parameters/classification/parameter_classification.h>
#include <chargefw/parameters/models/parameter_set.h>

#include <cstddef>
#include <optional>
#include <vector>

namespace chargefw::methods {

class Method;

struct ApplicableMethod {
    const Method* method = nullptr;
    const parameters::ParameterSet* parameter_set = nullptr;

    MethodOptions method_options;
    std::vector<parameters::ParameterClassification> classifications;

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
