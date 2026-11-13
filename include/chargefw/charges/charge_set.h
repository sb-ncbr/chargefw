#pragma once

#include <chargefw/charges/atomic_charges.h>

#include <cstddef>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace chargefw::charges {

struct ChargeTarget {
    std::size_t molecule_index = 0;
    std::optional<std::size_t> conformer_index;
};

struct ChargeAssignment {
    ChargeTarget target;
    AtomicCharges charges;
};

class ChargeSet {
  public:
    explicit ChargeSet(std::string method_id, std::vector<ChargeAssignment> assignments,
                       std::optional<std::string> parameter_set_id = std::nullopt);

    [[nodiscard]] auto method_id() const noexcept -> std::string_view;

    [[nodiscard]] auto parameter_set_id() const noexcept -> std::optional<std::string_view>;

    [[nodiscard]] auto assignments() const noexcept -> std::span<const ChargeAssignment>;

    [[nodiscard]] auto assignment(std::size_t index) const -> const ChargeAssignment&;

    [[nodiscard]] auto size() const noexcept -> std::size_t;
    [[nodiscard]] auto empty() const noexcept -> bool;

  private:
    std::string method_id_;
    std::optional<std::string> parameter_set_id_;
    std::vector<ChargeAssignment> assignments_;
};

} // namespace chargefw::charges
