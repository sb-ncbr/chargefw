#include <chargefw/core/bond.h>

#include <stdexcept>

#include <snitch/snitch.hpp>

namespace core = chargefw::core;

TEST_CASE("bond order converts to and from its numeric value", "[core][bond]") {
    CHECK(core::bond_order_value(core::BondOrder::SINGLE) == 1);
    CHECK(core::bond_order_value(core::BondOrder::DOUBLE) == 2);
    CHECK(core::bond_order_value(core::BondOrder::TRIPLE) == 3);
    CHECK(core::to_string(core::BondOrder::SINGLE) == "1");
    CHECK(core::to_string(core::BondOrder::DOUBLE) == "2");
    CHECK(core::to_string(core::BondOrder::TRIPLE) == "3");
    CHECK(core::bond_order_from_value(1) == core::BondOrder::SINGLE);
    CHECK(core::bond_order_from_value(2) == core::BondOrder::DOUBLE);
    CHECK(core::bond_order_from_value(3) == core::BondOrder::TRIPLE);
}

TEST_CASE("bond rejects invalid endpoints and orders", "[core][bond]") {
    CHECK_THROWS_AS((core::Bond{0, 0, core::BondOrder::SINGLE}), std::invalid_argument);
    CHECK_THROWS_AS(core::bond_order_from_value(4), std::invalid_argument);
    // Deliberately outside BondOrder, to check that it is rejected.
    // NOLINTNEXTLINE(clang-analyzer-optin.core.EnumCastOutOfRange)
    CHECK_THROWS_AS(core::to_string(static_cast<core::BondOrder>(4)), std::invalid_argument);
}
