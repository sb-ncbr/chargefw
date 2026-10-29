#pragma once

#include <chargefw/core/position.h>

#include <cstddef>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace chargefw::core {

class Conformer {
  public:
    explicit Conformer(std::vector<Position> positions, std::string name = {});

    [[nodiscard]] auto positions() const noexcept -> std::span<const Position>;

    [[nodiscard]] auto size() const noexcept -> std::size_t;
    [[nodiscard]] auto empty() const noexcept -> bool;

    [[nodiscard]] auto name() const noexcept -> std::string_view;

    [[nodiscard]] auto operator[](std::size_t index) const noexcept -> const Position&;
    [[nodiscard]] auto at(std::size_t index) const -> const Position&;

    [[nodiscard]] auto operator==(const Conformer&) const -> bool = default;

  private:
    std::vector<Position> positions_;
    std::string name_;
};

} // namespace chargefw::core
