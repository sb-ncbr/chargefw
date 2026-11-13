#include <chargefw/parameters/classification/classification_result.h>

#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace chargefw::parameters {

auto parameter_matching_from_string(const std::string_view value) -> ParameterMatching {
    if (value == "strict") {
        return ParameterMatching::strict;
    }
    if (value == "permissive") {
        return ParameterMatching::permissive;
    }
    throw std::invalid_argument{"parameter matching must be 'strict' or 'permissive', not '" +
                                std::string{value} + "'"};
}

auto to_string(const ParameterMatching value) -> std::string_view {
    switch (value) {
    case ParameterMatching::strict:
        return "strict";
    case ParameterMatching::permissive:
        return "permissive";
    }
    throw std::invalid_argument{"unknown parameter matching"};
}

ClassificationResult::ClassificationResult(ParameterClassification classification)
    : classification_{std::move(classification)} {}

ClassificationResult::ClassificationResult(std::vector<ClassificationIssue> issues)
    : issues_{std::move(issues)} {}

ClassificationResult::operator bool() const noexcept {
    return issues_.empty();
}

auto ClassificationResult::classification() const -> const ParameterClassification& {
    if (!issues_.empty()) {
        throw std::logic_error{"classification is not available for a failed result"};
    }

    return classification_;
}

auto ClassificationResult::issues() const noexcept -> const std::vector<ClassificationIssue>& {
    return issues_;
}

} // namespace chargefw::parameters
