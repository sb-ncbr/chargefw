#include <chargefw/calculation/calculation.h>
#include <chargefw/calculation/execution_policy.h>

#include <limits>
#include <optional>
#include <stdexcept>
#include <utility>

#include <snitch/snitch.hpp>

namespace calculation = chargefw::calculation;

TEST_CASE("execution selection parses from strings", "[calculation][execution-policy]") {
    CHECK(calculation::execution_selection_kind_from_string("auto") ==
          calculation::ExecutionSelectionKind::automatic);
    CHECK(calculation::execution_selection_kind_from_string("full") ==
          calculation::ExecutionSelectionKind::full);
    CHECK(calculation::execution_selection_kind_from_string("cutoff") ==
          calculation::ExecutionSelectionKind::cutoff);
    CHECK(calculation::execution_selection_kind_from_string("cover") ==
          calculation::ExecutionSelectionKind::cover);

    const auto bad_selection = [] {
        static_cast<void>(calculation::execution_selection_kind_from_string("unknown"));
    };
    CHECK_THROWS_AS(bad_selection(), std::invalid_argument);
}

TEST_CASE("execution policy validates mode and radius", "[calculation][execution-policy]") {
    using Mode = calculation::ExecutionMode;
    using Case = std::pair<Mode, std::optional<double>>;
    constexpr auto nan = std::numeric_limits<double>::quiet_NaN();
    constexpr auto infinity = std::numeric_limits<double>::infinity();

    const calculation::ExecutionPolicy default_policy;
    CHECK(default_policy.mode() == Mode::full);
    CHECK_FALSE(default_policy.radius().has_value());

    for (const auto& [mode, radius] :
         {Case{Mode::full, std::nullopt}, Case{Mode::cutoff, 8.0}, Case{Mode::cover, 12.0}}) {
        const calculation::ExecutionPolicy policy{mode, radius};
        CHECK(policy.mode() == mode);
        CHECK(policy.radius() == radius);
    }

    for (const auto& [mode, radius] :
         {Case{Mode::full, 8.0}, Case{Mode::cutoff, std::nullopt}, Case{Mode::cutoff, nan},
          Case{Mode::cutoff, infinity}, Case{Mode::cover, 7.99}}) {
        CAPTURE(calculation::to_string(mode), radius.value_or(-1.0));
        CHECK_THROWS_AS(static_cast<void>(calculation::ExecutionPolicy(mode, radius)),
                        std::invalid_argument);
    }
}

TEST_CASE("execution selection validates kind and radius", "[calculation][execution-policy]") {
    using Kind = calculation::ExecutionSelectionKind;
    using Case = std::pair<Kind, std::optional<double>>;

    const calculation::ExecutionSelection default_selection;
    CHECK(default_selection.kind() == Kind::automatic);
    CHECK_FALSE(default_selection.radius().has_value());

    for (const auto& [kind, radius] :
         {Case{Kind::automatic, 8.0}, Case{Kind::full, std::nullopt}}) {
        const calculation::ExecutionSelection selection{kind, radius};
        CHECK(selection.kind() == kind);
        CHECK(selection.radius() == radius);
    }

    for (const auto& [kind, radius] :
         {Case{Kind::automatic, 7.99}, Case{Kind::full, 8.0}, Case{Kind::cutoff, std::nullopt}}) {
        CAPTURE(calculation::to_string(kind), radius.value_or(-1.0));
        CHECK_THROWS_AS(static_cast<void>(calculation::ExecutionSelection(kind, radius)),
                        std::invalid_argument);
    }
}

TEST_CASE("resource policy exposes thresholds", "[calculation][execution-policy]") {
    const calculation::ResourcePolicy default_resources;
    CHECK(default_resources.cutoff_atom_threshold ==
          std::optional<std::size_t>{calculation::default_cutoff_atom_threshold});
    CHECK(default_resources.cover_atom_threshold ==
          std::optional<std::size_t>{calculation::default_cover_atom_threshold});

    const calculation::ResourcePolicy finite_resources{.cutoff_atom_threshold = 42,
                                                       .cover_atom_threshold = 84};
    CHECK(finite_resources.cutoff_atom_threshold == std::optional<std::size_t>{42});
    CHECK(finite_resources.cover_atom_threshold == std::optional<std::size_t>{84});

    const calculation::ResourcePolicy unlimited_resources{.cutoff_atom_threshold = std::nullopt,
                                                          .cover_atom_threshold = std::nullopt};
    CHECK_FALSE(unlimited_resources.cutoff_atom_threshold.has_value());
    CHECK_FALSE(unlimited_resources.cover_atom_threshold.has_value());

    CHECK(calculation::default_automatic_reduced_radius == 12.0);
}
