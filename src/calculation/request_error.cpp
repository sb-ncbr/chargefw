#include <chargefw/calculation/request_error.h>

#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace chargefw::calculation {

auto to_string(const RequestErrorKind value) -> std::string_view {
    switch (value) {
    case RequestErrorKind::unknown_method:
        return "unknown_method";
    case RequestErrorKind::unknown_parameter_set:
        return "unknown_parameter_set";
    case RequestErrorKind::incompatible_parameter_set:
        return "incompatible_parameter_set";
    case RequestErrorKind::invalid_execution:
        return "invalid_execution";
    }
    throw std::logic_error{"unknown request error kind"};
}

RequestError::RequestError(const RequestErrorKind kind, const std::string& reason,
                           std::optional<std::string> option)
    : std::invalid_argument{reason}, kind_{kind}, option_{std::move(option)} {}

auto RequestError::kind() const noexcept -> RequestErrorKind {
    return kind_;
}

auto RequestError::reason() const noexcept -> std::string_view {
    return what();
}

auto RequestError::option() const noexcept -> const std::optional<std::string>& {
    return option_;
}

} // namespace chargefw::calculation
