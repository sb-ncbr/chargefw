#pragma once

#include <cstddef>
#include <optional>
#include <vector>

namespace chargefw::core {

// Where in a molecule an issue is, in zero-based indices. The molecule itself is identified by
// the issue, because some issues always have one and others need not.
struct IssueLocation {
    std::optional<std::size_t> conformer_index = std::nullopt;
    // The atom or bond the issue's message names.
    std::optional<std::size_t> atom_index = std::nullopt;
    std::optional<std::size_t> bond_index = std::nullopt;
    // Every atom and bond the issue concerns, in index order, so that programs can find them
    // without reading the message.
    std::vector<std::size_t> atom_indices = {};
    std::vector<std::size_t> bond_indices = {};

    [[nodiscard]] auto operator==(const IssueLocation&) const -> bool = default;
};

} // namespace chargefw::core
