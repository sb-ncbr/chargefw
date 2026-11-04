#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

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
    // The conformer whose geometry the issue concerns, if it concerns one.
    std::optional<std::size_t> conformer_index = std::nullopt;
    // The atom the issue points to, if it concerns atoms: for atoms_too_close, the second atom of
    // the closest pair. The message names both.
    std::optional<std::size_t> atom_index = std::nullopt;
    // Every atom the issue concerns, in atom order: for atoms_too_close, each atom of a close pair.
    std::vector<std::size_t> atom_indices = {};
};

} // namespace chargefw::calculation
