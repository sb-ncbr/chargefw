#pragma once

#include <cstdint>
#include <string>
#include <string_view>

namespace chargefw::parameters {

enum class AtomParameterClassificationKind : std::uint8_t {
    PLAIN,
    HIGHEST_BOND_ORDER,
    BONDED_ELEMENTS
};

enum class BondParameterClassificationKind : std::uint8_t { PLAIN, BOND_ORDER };

struct AtomParameterKey {
    int atomic_number = 0;
    AtomParameterClassificationKind classification = AtomParameterClassificationKind::PLAIN;
    std::string type = "*";

    [[nodiscard]] auto operator==(const AtomParameterKey&) const -> bool = default;
};

struct BondTypeKey {
    BondParameterClassificationKind classification = BondParameterClassificationKind::PLAIN;
    std::string type = "*";

    [[nodiscard]] auto operator==(const BondTypeKey&) const -> bool = default;
};

// Equality compares the endpoints in their stored order.
struct BondParameterKey {
    AtomParameterKey first_atom;
    AtomParameterKey second_atom;
    BondTypeKey bond;

    [[nodiscard]] auto operator==(const BondParameterKey&) const -> bool = default;
};

[[nodiscard]] auto atom_classification_kind_from_string(std::string_view value)
    -> AtomParameterClassificationKind;

[[nodiscard]] auto bond_classification_kind_from_string(std::string_view value)
    -> BondParameterClassificationKind;

[[nodiscard]] auto to_string(AtomParameterClassificationKind kind) noexcept -> std::string_view;
[[nodiscard]] auto to_string(BondParameterClassificationKind kind) noexcept -> std::string_view;

} // namespace chargefw::parameters