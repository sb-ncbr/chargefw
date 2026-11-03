#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <string_view>

namespace chargefw::calculation {

enum class InputIssueKind : std::uint8_t {
    // The molecule contains carbon but no hydrogen atoms.
    hydrogens_missing,
};

[[nodiscard]] auto to_string(InputIssueKind value) -> std::string_view;

// A caution about an input molecule, such as hydrogens that were probably not added. Input issues
// do not exclude any method; assessment finds them once, independently of the plan.
struct InputIssue {
    InputIssueKind kind;
    std::string message;
    std::size_t molecule_index; // Zero-based in the supplied molecule collection.
    // The conformer whose geometry the issue concerns, if it concerns one.
    std::optional<std::size_t> conformer_index = std::nullopt;
};

} // namespace chargefw::calculation
