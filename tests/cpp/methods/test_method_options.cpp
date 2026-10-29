#include <chargefw/methods/method_options.h>

#include <array>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>

#include <snitch/snitch.hpp>

namespace methods = chargefw::methods;

TEST_CASE("method option types convert to stable strings", "[methods][method-options]") {
    CHECK(methods::to_string(methods::MethodOptionType::boolean) == "boolean");
    CHECK(methods::to_string(methods::MethodOptionType::floating_point) == "floating_point");
}

TEST_CASE("method options store defaults and support overrides", "[methods][method-options]") {
    const std::array schema{
        methods::MethodOptionSpec{.id = "enabled",
                                  .description = "Enable the method",
                                  .type = methods::MethodOptionType::boolean,
                                  .default_value = true,
                                  .choices = {}},
        methods::MethodOptionSpec{.id = "limit",
                                  .description = "Iteration limit",
                                  .type = methods::MethodOptionType::integer,
                                  .default_value = 25,
                                  .choices = {}},
        methods::MethodOptionSpec{.id = "radius",
                                  .description = "Search radius",
                                  .type = methods::MethodOptionType::floating_point,
                                  .default_value = 1.5,
                                  .choices = {},
                                  .minimum = 0.0,
                                  .minimum_inclusive = false},
        methods::MethodOptionSpec{.id = "label",
                                  .description = "Label",
                                  .type = methods::MethodOptionType::string,
                                  .default_value = std::string{"default"},
                                  .choices = {}}};

    auto options = methods::make_default_options(schema);

    CHECK(options.contains("enabled"));
    CHECK(options.contains("limit"));
    CHECK(options.contains("radius"));
    CHECK(options.contains("label"));
    CHECK(options.get<bool>("enabled"));
    CHECK(options.get<int>("limit") == 25);
    CHECK(options.get<double>("radius") == 1.5);
    CHECK(options.get<std::string>("label") == "default");

    options.set("limit", 50);
    CHECK(options.get<int>("limit") == 50);

    CHECK_NOTHROW(methods::validate_method_options(schema, options));
}

TEST_CASE("method options reject missing, wrong type, and out-of-bounds values",
          "[methods][method-options]") {
    const std::array schema{
        methods::MethodOptionSpec{.id = "enabled",
                                  .description = "Enable the method",
                                  .type = methods::MethodOptionType::boolean,
                                  .default_value = true,
                                  .choices = {}},
        methods::MethodOptionSpec{.id = "limit",
                                  .description = "Iteration limit",
                                  .type = methods::MethodOptionType::integer,
                                  .default_value = 25,
                                  .choices = {},
                                  .maximum = 100},
        methods::MethodOptionSpec{.id = "radius",
                                  .description = "Search radius",
                                  .type = methods::MethodOptionType::floating_point,
                                  .default_value = 1.5,
                                  .choices = {},
                                  .minimum = 0.0,
                                  .minimum_inclusive = false},
        methods::MethodOptionSpec{.id = "label",
                                  .description = "Label",
                                  .type = methods::MethodOptionType::string,
                                  .default_value = std::string{"default"},
                                  .choices = {}}};

    auto options = methods::make_default_options(schema);

    // missing required option
    CHECK_THROWS_AS(methods::validate_method_options(schema, methods::MethodOptions{}),
                    std::invalid_argument);

    // wrong type for option
    const auto wrong_type = [&] {
        auto invalid = methods::make_default_options(schema);
        invalid.set("limit", 50.0);
        methods::validate_method_options(schema, invalid);
    };
    CHECK_THROWS_AS(wrong_type(), std::invalid_argument);

    // value below exclusive minimum
    const auto below_minimum = [&] {
        auto invalid = methods::make_default_options(schema);
        invalid.set("radius", 0.0);
        methods::validate_method_options(schema, invalid);
    };
    CHECK_THROWS_AS(below_minimum(), std::invalid_argument);

    // nonfinite values, including infinity above a lone lower bound
    for (const auto value :
         {std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity(),
          -std::numeric_limits<double>::infinity()}) {
        CAPTURE(value);
        auto invalid = methods::make_default_options(schema);
        invalid.set("radius", value);
        CHECK_THROWS_AS(methods::validate_method_options(schema, invalid), std::invalid_argument);
    }

    // inclusive maximum accepts the bound and rejects values above it
    auto at_maximum = methods::make_default_options(schema);
    at_maximum.set("limit", 100);
    CHECK_NOTHROW(methods::validate_method_options(schema, at_maximum));
    const auto above_maximum = [&] {
        auto invalid = methods::make_default_options(schema);
        invalid.set("limit", 101);
        methods::validate_method_options(schema, invalid);
    };
    CHECK_THROWS_AS(above_maximum(), std::invalid_argument);

    // unknown option key
    const auto unknown_key = [&] {
        auto invalid = methods::make_default_options(schema);
        invalid.set("extra", true);
        methods::validate_method_options(schema, invalid);
    };
    CHECK_THROWS_AS(unknown_key(), std::invalid_argument);

    // unknown option lookup throws out_of_range
    CHECK_THROWS_AS(options.get<int>("missing"), std::out_of_range);

    // wrong type lookup throws invalid_argument
    CHECK_THROWS_AS(options.get<double>("limit"), std::invalid_argument);
}

TEST_CASE("method option schema rejects mismatched bounds and duplicate IDs",
          "[methods][method-options]") {
    const std::array mismatched_bound_schema{
        methods::MethodOptionSpec{.id = "limit",
                                  .description = "Iteration limit",
                                  .type = methods::MethodOptionType::integer,
                                  .default_value = 1,
                                  .choices = {},
                                  .minimum = 0.0}};

    CHECK_THROWS_AS(methods::validate_method_option_schema(mismatched_bound_schema),
                    std::invalid_argument);

    const std::array duplicate_id_schema{
        methods::MethodOptionSpec{.id = "limit",
                                  .description = "First limit",
                                  .type = methods::MethodOptionType::integer,
                                  .default_value = 25,
                                  .choices = {}},
        methods::MethodOptionSpec{.id = "limit",
                                  .description = "Second limit",
                                  .type = methods::MethodOptionType::integer,
                                  .default_value = 50,
                                  .choices = {}}};

    CHECK_THROWS_AS(methods::validate_method_option_schema(duplicate_id_schema),
                    std::invalid_argument);
}

TEST_CASE("method option schema rejects wrong default and choice types",
          "[methods][method-options]") {
    const std::array wrong_default_type_schema{
        methods::MethodOptionSpec{.id = "limit",
                                  .description = "Iteration limit",
                                  .type = methods::MethodOptionType::integer,
                                  .default_value = 25.0,
                                  .choices = {}}};

    CHECK_THROWS_AS(methods::validate_method_option_schema(wrong_default_type_schema),
                    std::invalid_argument);

    const std::array wrong_choice_type_schema{
        methods::MethodOptionSpec{.id = "limit",
                                  .description = "Iteration limit",
                                  .type = methods::MethodOptionType::integer,
                                  .default_value = 25,
                                  .choices = {25, 50.0}}};

    CHECK_THROWS_AS(methods::validate_method_option_schema(wrong_choice_type_schema),
                    std::invalid_argument);

    const std::array default_not_allowed_schema{
        methods::MethodOptionSpec{.id = "limit",
                                  .description = "Iteration limit",
                                  .type = methods::MethodOptionType::integer,
                                  .default_value = 25,
                                  .choices = {50, 100}}};

    CHECK_THROWS_AS(methods::validate_method_option_schema(default_not_allowed_schema),
                    std::invalid_argument);
}

TEST_CASE("method option choices restrict allowed values", "[methods][method-options]") {
    const std::array choice_schema{
        methods::MethodOptionSpec{.id = "mode",
                                  .description = "Calculation mode",
                                  .type = methods::MethodOptionType::string,
                                  .default_value = std::string{"fast"},
                                  .choices = {std::string{"fast"}, std::string{"accurate"}}},
        methods::MethodOptionSpec{.id = "limit",
                                  .description = "Iteration limit",
                                  .type = methods::MethodOptionType::integer,
                                  .default_value = 25,
                                  .choices = {10, 25, 50}}};

    auto choice_options = methods::make_default_options(choice_schema);
    CHECK_NOTHROW(methods::validate_method_options(choice_schema, choice_options));

    // string choice not in allowed set
    const auto bad_string_choice = [&] {
        auto invalid = methods::make_default_options(choice_schema);
        invalid.set("mode", std::string{"slow"});
        methods::validate_method_options(choice_schema, invalid);
    };
    CHECK_THROWS_AS(bad_string_choice(), std::invalid_argument);

    // integer choice not in allowed set
    const auto bad_int_choice = [&] {
        auto invalid = methods::make_default_options(choice_schema);
        invalid.set("limit", 100);
        methods::validate_method_options(choice_schema, invalid);
    };
    CHECK_THROWS_AS(bad_int_choice(), std::invalid_argument);
}
