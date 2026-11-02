#pragma once

#include <compare>

namespace chargefw::core {

struct Position {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    [[nodiscard]] auto operator==(const Position&) const -> bool = default;
    // Lexicographic floating-point comparison by x, then y, then z; not an order by distance. Like
    // ==, it treats -0.0 and 0.0 as equal. It may be unordered when a compared coordinate is NaN.
    [[nodiscard]] auto operator<=>(const Position&) const = default;
};

[[nodiscard]] auto squared_distance(const Position& first, const Position& second) noexcept
    -> double;

[[nodiscard]] auto distance(const Position& first, const Position& second) noexcept -> double;

} // namespace chargefw::core
