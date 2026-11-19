#include <chargefw/io/parse_error.h>

#include <utility>

namespace chargefw::io {
namespace {

[[nodiscard]] auto describe(const std::string& reason, const ParseLocation& location)
    -> std::string {
    auto text = std::string{};
    if (location.line.has_value()) {
        text = (location.source_name.empty() ? "line " : location.source_name + ":") +
               std::to_string(*location.line) + ": ";
    } else if (!location.source_name.empty()) {
        text = location.source_name + ": ";
    }
    if (location.record_index.has_value()) {
        text +=
            "record " + std::to_string(*location.record_index + 1) +
            (location.record_id.empty() ? "" : " ('" + location.record_id.display_string() + "')") +
            ": ";
    }
    return text + reason;
}

} // namespace

ParseError::ParseError(std::string reason, ParseLocation location)
    : std::runtime_error{describe(reason, location)}, reason_{std::move(reason)},
      location_{std::move(location)} {}

auto ParseError::reason() const noexcept -> const std::string& {
    return reason_;
}

auto ParseError::location() const noexcept -> const ParseLocation& {
    return location_;
}

} // namespace chargefw::io
