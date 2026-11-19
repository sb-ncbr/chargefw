#pragma once

#include <cstddef>
#include <span>
#include <stdexcept>
#include <vector>

namespace chargefw::charges {

// Thrown by AtomicCharges for values that are not finite numbers.
class NonFiniteCharges : public std::invalid_argument {
  public:
    using std::invalid_argument::invalid_argument;
};

class AtomicCharges {
  public:
    explicit AtomicCharges(std::vector<double> values);

    [[nodiscard]] auto values() const noexcept -> std::span<const double>;

    [[nodiscard]] auto size() const noexcept -> std::size_t;
    [[nodiscard]] auto empty() const noexcept -> bool;

    [[nodiscard]] auto total() const noexcept -> double;

    [[nodiscard]] auto operator[](std::size_t index) const noexcept -> double;
    [[nodiscard]] auto at(std::size_t index) const -> double;

  private:
    std::vector<double> values_;
};

} // namespace chargefw::charges