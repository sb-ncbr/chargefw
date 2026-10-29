#pragma once

namespace chargefw::core {

struct Position {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    [[nodiscard]] auto operator==(const Position&) const -> bool = default;
};

[[nodiscard]] auto squared_distance(const Position& first, const Position& second) noexcept
    -> double;

[[nodiscard]] auto distance(const Position& first, const Position& second) noexcept -> double;

} // namespace chargefw::core
