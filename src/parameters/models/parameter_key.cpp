#include <chargefw/parameters/models/parameter_key.h>

#include <stdexcept>
#include <string>

namespace chargefw::parameters {

auto atom_classification_kind_from_string(const std::string_view value)
    -> AtomParameterClassificationKind {
    if (value == "plain") {
        return AtomParameterClassificationKind::plain;
    }

    if (value == "hbo") {
        return AtomParameterClassificationKind::highest_bond_order;
    }

    if (value == "bonded") {
        return AtomParameterClassificationKind::bonded_elements;
    }

    throw std::invalid_argument{"unknown atom parameter classification '" + std::string{value} +
                                "'"};
}

auto bond_classification_kind_from_string(const std::string_view value)
    -> BondParameterClassificationKind {
    if (value == "plain") {
        return BondParameterClassificationKind::plain;
    }

    if (value == "bo") {
        return BondParameterClassificationKind::bond_order;
    }

    throw std::invalid_argument{"unknown bond parameter classification '" + std::string{value} +
                                "'"};
}

auto to_string(const AtomParameterClassificationKind kind) noexcept -> std::string_view {
    switch (kind) {
    case AtomParameterClassificationKind::plain:
        return "plain";

    case AtomParameterClassificationKind::highest_bond_order:
        return "hbo";

    case AtomParameterClassificationKind::bonded_elements:
        return "bonded";
    }

    return "unknown";
}

auto to_string(const BondParameterClassificationKind kind) noexcept -> std::string_view {
    switch (kind) {
    case BondParameterClassificationKind::plain:
        return "plain";

    case BondParameterClassificationKind::bond_order:
        return "bo";
    }

    return "unknown";
}

} // namespace chargefw::parameters