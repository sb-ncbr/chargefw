#include <chargefw/methods/method_failure.h>

namespace chargefw::methods {

auto to_string(const FailureKind value) -> std::string_view {
    switch (value) {
    case FailureKind::no_solution:
        return "no_solution";
    case FailureKind::not_converged:
        return "not_converged";
    case FailureKind::invalid_parameter_value:
        return "invalid_parameter_value";
    }
    throw std::invalid_argument{"unknown method failure kind"};
}

MethodFailure::MethodFailure(const FailureKind kind, const std::string& message)
    : std::runtime_error{message}, kind_{kind} {}

auto MethodFailure::kind() const noexcept -> FailureKind {
    return kind_;
}

} // namespace chargefw::methods
