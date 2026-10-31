#pragma once

#include <cstddef>
#include <optional>
#include <span>
#include <string>
#include <string_view>

#include <chargefw/adapters/molecule_record.h>

namespace chargefw::adapters::detail {

// A source label, or nothing when it is absent, empty, or an mmCIF placeholder ('.' or '?').
[[nodiscard]] inline auto known_label(const std::optional<std::string>& value)
    -> std::optional<std::string_view> {
    if (!value || value->empty() || *value == "." || *value == "?") {
        return std::nullopt;
    }
    return *value;
}

// Chemical identity prefers label names, which follow the Chemical Component Dictionary; locations
// shown to a person prefer author names, which PDB files and viewers use.
[[nodiscard]] inline auto first_known_label(const std::optional<std::string>& preferred,
                                            const std::optional<std::string>& alternate,
                                            const std::string_view fallback) -> std::string_view {
    return known_label(preferred)
        .or_else([&alternate] { return known_label(alternate); })
        .value_or(fallback);
}

// A residue as users know it from PDB files, such as "MSE A 52A".
[[nodiscard]] inline auto
residue_description(const std::string_view name, const std::string_view chain,
                    const std::string_view sequence, const std::string_view insertion_code)
    -> std::string {
    auto result = std::string{name} + ' ' + std::string{chain} + ' ' + std::string{sequence};
    if (!insertion_code.empty() && insertion_code != "." && insertion_code != "?") {
        result += insertion_code;
    }
    return result;
}

struct ComponentAtom {
    std::string name;
    std::size_t index;
};

[[nodiscard]] inline auto
canonical_component_name(const std::optional<SourceStructuralLabels>& labels,
                         const std::string_view fallback) -> std::string_view {
    return labels ? first_known_label(labels->label.residue, labels->author.residue, fallback)
                  : fallback;
}

[[nodiscard]] inline auto canonical_atom_name(const std::optional<SourceStructuralLabels>& labels,
                                              const std::string_view fallback) -> std::string_view {
    return labels ? first_known_label(labels->label.atom, labels->author.atom, fallback) : fallback;
}

struct ComponentView {
    std::string_view name;
    std::span<const ComponentAtom> atoms;

    [[nodiscard]] auto find_atom(const std::string_view atom_name) const
        -> std::optional<std::size_t> {
        for (const auto& atom : atoms) {
            if (atom.name == atom_name) {
                return atom.index;
            }
        }
        return std::nullopt;
    }
};

} // namespace chargefw::adapters::detail
