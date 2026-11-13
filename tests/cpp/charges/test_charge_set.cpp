#include <chargefw/charges/charge_set.h>

#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>

#include <snitch/snitch.hpp>

namespace charges = chargefw::charges;

TEST_CASE("charge set stores method, parameter, and assignments", "[charges][charge-set]") {
    const charges::ChargeTarget conformer_target{.molecule_index = 0, .conformer_index = 0};
    const charges::ChargeAssignment water_assignment{
        .target = conformer_target, .charges = charges::AtomicCharges{{-0.8, 0.4, 0.4}}};
    const charges::ChargeAssignment pair_assignment{
        .target = charges::ChargeTarget{.molecule_index = 1, .conformer_index = std::nullopt},
        .charges = charges::AtomicCharges{{1.0, -1.0}}};

    const charges::ChargeSet charge_set{
        "formal", {water_assignment, pair_assignment}, std::string{"default"}};

    CHECK(charge_set.method_id() == std::string_view{"formal"});
    REQUIRE(charge_set.parameter_set_id().has_value());
    CHECK(*charge_set.parameter_set_id() == std::string_view{"default"});
    REQUIRE(charge_set.size() == 2);
    CHECK_FALSE(charge_set.empty());
    REQUIRE(charge_set.assignments().size() == 2);
    CHECK(charge_set.assignment(0).target.molecule_index == 0);
    CHECK(charge_set.assignment(1).target.molecule_index == 1);
}

TEST_CASE("parameterless charge set has no parameter-set id", "[charges][charge-set]") {
    const charges::ChargeSet parameterless_charge_set{"dummy", {}};

    CHECK_FALSE(parameterless_charge_set.parameter_set_id().has_value());
    CHECK(parameterless_charge_set.empty());
}

TEST_CASE("charge set rejects empty method id and invalid assignment index",
          "[charges][charge-set]") {
    const charges::ChargeTarget conformer_target{.molecule_index = 0, .conformer_index = 0};
    const charges::ChargeAssignment water_assignment{
        .target = conformer_target, .charges = charges::AtomicCharges{{-0.8, 0.4, 0.4}}};
    const charges::ChargeAssignment pair_assignment{
        .target = charges::ChargeTarget{.molecule_index = 1, .conformer_index = std::nullopt},
        .charges = charges::AtomicCharges{{1.0, -1.0}}};
    const charges::ChargeSet charge_set{
        "formal", {water_assignment, pair_assignment}, std::string{"default"}};
    CHECK_THROWS_AS((charges::ChargeSet{"", {}}), std::invalid_argument);
    CHECK_THROWS_AS(charge_set.assignment(2), std::out_of_range);
}
