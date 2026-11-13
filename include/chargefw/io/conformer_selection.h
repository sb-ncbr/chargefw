#pragma once

#include <cstdint>
#include <string_view>

namespace chargefw::io {

enum class ConformerSelection : std::uint8_t {
    first,
    all,
};

[[nodiscard]] auto conformer_selection_from_string(std::string_view value) -> ConformerSelection;
[[nodiscard]] auto to_string(ConformerSelection value) -> std::string_view;

} // namespace chargefw::io
