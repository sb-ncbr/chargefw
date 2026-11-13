#include "calculation/calculation_request.h"
#include "methods/method_applicability.h"
#include "support/test_methods.h"
#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <chargefw/calculation/calculation.h>
#include <chargefw/charges/charge_set.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/calculation_input.h>
#include <chargefw/methods/method.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/methods/method_registry.h>
#include <chargefw/parameters/classification/parameter_classification.h>
#include <chargefw/parameters/models/atom_parameters.h>
#include <chargefw/parameters/models/parameter_key.h>
#include <chargefw/parameters/models/parameter_set.h>
#include <chargefw/parameters/models/parameter_set_metadata.h>

#include <cstddef>
#include <optional>
#include <snitch/snitch.hpp>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace calculation = chargefw::calculation;
namespace charges = chargefw::charges;
namespace core = chargefw::core;
namespace features = chargefw::features;
namespace methods = chargefw::methods;
namespace parameters = chargefw::parameters;

namespace {

class WrongSizeMethod final : public chargefw::test::StubMethod {
  public:
    WrongSizeMethod() : StubMethod{"wrong-size-test"} {}

    [[nodiscard]] auto calculate(const methods::CalculationInput& /* unused */) const
        -> charges::AtomicCharges override {
        return charges::AtomicCharges{std::vector{0.0}};
    }
};

class AtomParameterMethod final : public chargefw::test::StubMethod {
  public:
    AtomParameterMethod() : StubMethod{"atom-parameter-test", {.atom_parameters = {"value"}}} {}

    [[nodiscard]] auto calculate(const methods::CalculationInput& input) const
        -> charges::AtomicCharges override {
        const auto& molecule = input.molecule();
        const auto value = input.parameters().atom("value");

        std::vector<double> values;
        values.reserve(molecule.atom_count());

        for (std::size_t atom_index = 0; atom_index < molecule.atom_count(); ++atom_index) {
            values.push_back(value[atom_index]);
        }

        return charges::AtomicCharges{std::move(values)};
    }
};

class GeometryMethod final : public chargefw::test::StubMethod {
  public:
    GeometryMethod() : StubMethod{"geometry-test", {.coordinates = true}} {}

    [[nodiscard]] auto calculate(const methods::CalculationInput& input) const
        -> charges::AtomicCharges override {
        const auto& molecule = input.molecule();
        const auto& geometry = input.geometry();

        std::vector<double> values;
        values.reserve(molecule.atom_count());

        for (std::size_t atom_index = 0; atom_index < molecule.atom_count(); ++atom_index) {
            values.push_back(geometry.position(atom_index).x);
        }

        return charges::AtomicCharges{std::move(values)};
    }
};

auto make_collection() -> core::MoleculeCollection {
    std::vector molecules{chargefw::test::make_water(),
                          chargefw::test::make_formally_charged_pair()};

    return core::MoleculeCollection{std::move(molecules), "test-collection"};
}

auto make_parameter_set() -> parameters::ParameterSet {
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{
            .id = "test-parameters", .method_id = "atom-parameter-test", .name = "Test parameters"},
        {},
        parameters::AtomParameters{
            {{.key = chargefw::test::atom_key(
                  1, parameters::AtomParameterClassificationKind::bonded_elements, "O"),
              .parameters = {{.name = "value", .value = 1.0}}},
             {.key = chargefw::test::atom_key(
                  8, parameters::AtomParameterClassificationKind::bonded_elements, "HH"),
              .parameters = {{.name = "value", .value = 2.0}}},
             {.key = chargefw::test::atom_key(7, parameters::AtomParameterClassificationKind::plain,
                                              "*"),
              .parameters = {{.name = "value", .value = 3.0}}},
             {.key = chargefw::test::atom_key(
                  17, parameters::AtomParameterClassificationKind::plain, "*"),
              .parameters = {{.name = "value", .value = 4.0}}}}}};
}

auto make_parameterized_candidate(const AtomParameterMethod& method,
                                  const parameters::ParameterSet& parameter_set)
    -> methods::ApplicableMethod {
    return methods::ApplicableMethod{
        .method = &method,
        .parameter_set = &parameter_set,
        .method_options = methods::MethodOptions{},
        .classifications = {
            parameters::ParameterClassification{
                parameters::AtomParameterClassification{std::vector<std::size_t>{1, 0, 0}}},
            parameters::ParameterClassification{
                parameters::AtomParameterClassification{std::vector<std::size_t>{2, 3}}}}};
}

} // namespace

TEST_CASE("method calculation dispatches parameterized and geometry-dependent targets",
          "[methods][method-calculation]") {
    const auto collection = make_collection();
    const features::PreparedMoleculeCollection prepared{collection};

    const auto& registry = methods::method_registry();
    const auto* dummy = registry.find("dummy");

    REQUIRE(dummy != nullptr);

    const methods::ApplicableMethod dummy_candidate{
        .method = dummy,
        .parameter_set = nullptr,
        .method_options = methods::make_default_options(dummy->option_schema()),
        .classifications = {}};

    const auto dummy_charges =
        calculation::calculate({.molecules = prepared, .selected = dummy_candidate}).charges;

    CHECK(dummy_charges.method_id() == std::string_view{"dummy"});
    CHECK(!dummy_charges.parameter_set_id().has_value());
    REQUIRE(dummy_charges.size() == prepared.size());

    CHECK(dummy_charges.assignment(0).target.molecule_index == 0);
    CHECK(!dummy_charges.assignment(0).target.conformer_index.has_value());
    CHECK(dummy_charges.assignment(0).charges.size() == collection[0].atom_count());

    for (const auto value : dummy_charges.assignment(0).charges.values()) {
        CHECK(value == 0.0);
    }

    const AtomParameterMethod parameterized_method;
    const auto parameter_set = make_parameter_set();

    const auto parameterized_candidate =
        make_parameterized_candidate(parameterized_method, parameter_set);

    const auto parameterized_charges =
        calculation::calculate({.molecules = prepared, .selected = parameterized_candidate})
            .charges;

    CHECK(parameterized_charges.method_id() == std::string_view{"atom-parameter-test"});
    REQUIRE(parameterized_charges.parameter_set_id().has_value());
    CHECK(*parameterized_charges.parameter_set_id() == std::string_view{"test-parameters"});
    REQUIRE(parameterized_charges.size() == prepared.size());

    const auto& water_charges = parameterized_charges.assignment(0).charges;

    REQUIRE(water_charges.size() == 3);
    CHECK(water_charges[0] == 2.0);
    CHECK(water_charges[1] == 1.0);
    CHECK(water_charges[2] == 1.0);

    const auto& charged_pair_charges = parameterized_charges.assignment(1).charges;

    REQUIRE(charged_pair_charges.size() == 2);
    CHECK(charged_pair_charges[0] == 3.0);
    CHECK(charged_pair_charges[1] == 4.0);

    const auto null_method_call = [&] {
        const methods::ApplicableMethod invalid{.method = nullptr,
                                                .parameter_set = nullptr,
                                                .method_options = {},
                                                .classifications = {}};
        static_cast<void>(calculation::calculate({.molecules = prepared, .selected = invalid}));
    };
    CHECK_THROWS_AS(null_method_call(), std::invalid_argument);

    const auto missing_classification_call = [&] {
        const methods::ApplicableMethod invalid{.method = &parameterized_method,
                                                .parameter_set = &parameter_set,
                                                .method_options = {},
                                                .classifications = {}};
        static_cast<void>(calculation::calculate({.molecules = prepared, .selected = invalid}));
    };
    CHECK_THROWS_AS(missing_classification_call(), std::invalid_argument);

    const auto wrong_charge_count_call = [&] {
        const WrongSizeMethod wrong_size_method;
        const methods::ApplicableMethod invalid{.method = &wrong_size_method,
                                                .parameter_set = nullptr,
                                                .method_options = {},
                                                .classifications = {}};
        static_cast<void>(calculation::calculate({.molecules = prepared, .selected = invalid}));
    };
    CHECK_THROWS_AS(wrong_charge_count_call(), std::runtime_error);

    const auto invalid_classification_call = [&] {
        const methods::ApplicableMethod invalid{
            .method = &parameterized_method,
            .parameter_set = &parameter_set,
            .method_options = {},
            .classifications = {
                parameters::ParameterClassification{
                    parameters::AtomParameterClassification{std::vector<std::size_t>{999, 0, 0}}},
                parameters::ParameterClassification{
                    parameters::AtomParameterClassification{std::vector<std::size_t>{2, 3}}}}};
        static_cast<void>(calculation::calculate({.molecules = prepared, .selected = invalid}));
    };
    CHECK_THROWS_AS(invalid_classification_call(), std::invalid_argument);

    const GeometryMethod geometry_method;

    const auto two_conformer_collection = core::MoleculeCollection{
        std::vector{chargefw::test::make_two_conformer_water()}, "two-conformer-test"};

    const features::PreparedMoleculeCollection prepared_two_conformer_collection{
        two_conformer_collection};

    const methods::ApplicableMethod geometry_candidate{.method = &geometry_method,
                                                       .parameter_set = nullptr,
                                                       .method_options = {},
                                                       .classifications = {}};

    const auto geometry_charges =
        calculation::calculate(
            {.molecules = prepared_two_conformer_collection, .selected = geometry_candidate})
            .charges;

    CHECK(geometry_charges.method_id() == std::string_view{"geometry-test"});
    CHECK(!geometry_charges.parameter_set_id().has_value());
    REQUIRE(geometry_charges.size() == 2);

    CHECK(geometry_charges.assignment(0).target.molecule_index == 0);
    CHECK(geometry_charges.assignment(0).target.conformer_index == std::optional<std::size_t>{0});

    CHECK(geometry_charges.assignment(1).target.molecule_index == 0);
    CHECK(geometry_charges.assignment(1).target.conformer_index == std::optional<std::size_t>{1});

    REQUIRE(geometry_charges.assignment(0).charges.size() == 3);
    REQUIRE(geometry_charges.assignment(1).charges.size() == 3);

    CHECK(geometry_charges.assignment(0).charges[1] == 0.7570);
    CHECK(geometry_charges.assignment(1).charges[1] == 1.1000);
}
