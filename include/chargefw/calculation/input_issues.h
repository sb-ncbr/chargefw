#pragma once

#include <chargefw/core/issue_location.h>

#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>

namespace chargefw::calculation {

enum class InputIssueKind : std::uint8_t {
    // The molecule contains carbon but no hydrogen atoms.
    hydrogens_missing,
    // Atoms of a conformer are closer than 0.7 times the sum of their covalent radii. Reported
    // only when the executed method uses coordinates.
    atoms_too_close,
};

[[nodiscard]] auto to_string(InputIssueKind value) -> std::string_view;

// A caution about an input molecule, such as hydrogens that were probably not added. Input issues
// do not exclude any method. Assessment finds those independent of the plan; execution adds those
// that concern only the executed method.
struct InputIssue {
    InputIssueKind kind;
    std::string message;
    std::size_t molecule_index; // Zero-based in the supplied molecule collection.
    // For atoms_too_close: the conformer, the second atom of the closest pair (the message names
    // both), and every atom of a close pair.
    core::IssueLocation location = {};
};

} // namespace chargefw::calculation
