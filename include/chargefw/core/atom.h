#pragma once

#include <string>
#include <string_view>

namespace chargefw::core {
class Atom {
  public:
    explicit Atom(int atomic_number, int formal_charge = 0, std::string source_name = {});

    [[nodiscard]] auto atomic_number() const noexcept -> int;
    [[nodiscard]] auto formal_charge() const noexcept -> int;
    [[nodiscard]] auto name() const noexcept -> std::string_view;

    [[nodiscard]] auto operator==(const Atom&) const -> bool = default;

  private:
    int atomic_number_;
    int formal_charge_;
    std::string name_;
};
} // namespace chargefw::core