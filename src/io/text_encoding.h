#pragma once

#include <algorithm>
#include <cstddef>
#include <optional>
#include <string_view>

namespace chargefw::io::detail {

// Text formats other than JSON are read as ASCII. Returns the one-based line of the first byte
// outside ASCII, or nothing when the text is ASCII.
[[nodiscard]] inline auto first_non_ascii_line(const std::string_view text)
    -> std::optional<std::size_t> {
    const auto byte =
        std::ranges::find_if(text, [](const unsigned char value) { return value >= 0x80; });
    if (byte == text.end()) {
        return std::nullopt;
    }
    return 1 + static_cast<std::size_t>(std::ranges::count(text.begin(), byte, '\n'));
}

} // namespace chargefw::io::detail
