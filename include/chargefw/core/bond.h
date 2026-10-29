#pragma once

#include <cstddef>
#include <cstdint>
#include <string_view>

namespace chargefw::core {
enum class BondOrder : std::uint8_t { SINGLE = 1, DOUBLE = 2, TRIPLE = 3 };

[[nodiscard]] constexpr auto bond_order_value(const BondOrder order) noexcept -> std::uint8_t {
    return static_cast<std::uint8_t>(order);
}

[[nodiscard]] auto bond_order_from_value(int value) -> BondOrder;
[[nodiscard]] auto to_string(BondOrder order) -> std::string_view;

class Bond {
  public:
    explicit Bond(std::size_t first_atom_index, std::size_t second_atom_index,
                  BondOrder order = BondOrder::SINGLE);

    [[nodiscard]] auto first_atom_index() const noexcept -> std::size_t;
    [[nodiscard]] auto second_atom_index() const noexcept -> std::size_t;
    [[nodiscard]] auto order() const noexcept -> BondOrder;

    [[nodiscard]] auto operator==(const Bond&) const -> bool = default;

  private:
    std::size_t first_atom_index_;
    std::size_t second_atom_index_;
    BondOrder order_;
};

} // namespace chargefw::core
