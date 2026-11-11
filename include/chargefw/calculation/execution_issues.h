#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <string_view>

namespace chargefw::calculation {

enum class ExecutionIssueKind : std::uint8_t {
    resource_threshold_exceeded,
    unsupported_execution_mode,
};

[[nodiscard]] auto to_string(ExecutionIssueKind value) -> std::string_view;

// Why an execution mode of a method cannot be used, or a caution about using it: a warning of a
// plan, or a reason for a rejection.
struct ExecutionIssue {
    ExecutionIssueKind kind;
    std::string message;
    std::optional<std::size_t> molecule_index = std::nullopt;
};

} // namespace chargefw::calculation
