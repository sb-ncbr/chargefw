#include <chargefw/io/output_error.h>

#include <cstddef>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace chargefw::io {
namespace {

// Numbers count from one here and from zero in the fields.
[[nodiscard]] auto describe(const std::string& reason,
                            const std::optional<std::size_t> molecule_index,
                            const core::IssueLocation& location) -> std::string {
    auto text = std::string{};
    const auto add = [&text](const std::string_view name, const std::optional<std::size_t> index) {
        if (index.has_value()) {
            text +=
                (text.empty() ? "" : ", ") + std::string{name} + " " + std::to_string(*index + 1);
        }
    };
    add("molecule", molecule_index);
    add("conformer", location.conformer_index);
    if (location.atom_index.has_value() && location.other_atom_index.has_value()) {
        text += (text.empty() ? "" : ", ") + std::string{"atoms "} +
                std::to_string(*location.atom_index + 1) + " and " +
                std::to_string(*location.other_atom_index + 1);
    } else {
        add("atom", location.atom_index);
    }
    add("bond", location.bond_index);
    return text.empty() ? reason : text + ": " + reason;
}

} // namespace

auto to_string(const OutputErrorKind value) -> std::string_view {
    switch (value) {
    case OutputErrorKind::missing_coordinates:
        return "missing_coordinates";
    case OutputErrorKind::unsupported_shape:
        return "unsupported_shape";
    case OutputErrorKind::unrepresentable:
        return "unrepresentable";
    case OutputErrorKind::target_mismatch:
        return "target_mismatch";
    case OutputErrorKind::existing_charges:
        return "existing_charges";
    case OutputErrorKind::no_charges:
        return "no_charges";
    case OutputErrorKind::unsupported_source:
        return "unsupported_source";
    }
    throw std::logic_error{"unknown output error kind"};
}

OutputError::OutputError(const OutputErrorKind kind, std::string reason,
                         const std::optional<std::size_t> molecule_index,
                         core::IssueLocation location)
    : std::invalid_argument{describe(reason, molecule_index, location)}, kind_{kind},
      reason_{std::move(reason)}, molecule_index_{molecule_index}, location_{std::move(location)} {}

auto OutputError::kind() const noexcept -> OutputErrorKind {
    return kind_;
}

auto OutputError::reason() const noexcept -> const std::string& {
    return reason_;
}

auto OutputError::molecule_index() const noexcept -> std::optional<std::size_t> {
    return molecule_index_;
}

auto OutputError::location() const noexcept -> const core::IssueLocation& {
    return location_;
}

} // namespace chargefw::io
