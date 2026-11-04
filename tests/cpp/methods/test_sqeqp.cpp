#include "support/test_calculation.h"
#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/position.h>
#include <chargefw/features/conformer_features.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/methods/calculation_input.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/methods/method_registry.h>
#include <chargefw/parameters/classification/parameter_classification.h>
#include <chargefw/parameters/models/atom_parameters.h>
#include <chargefw/parameters/models/bond_parameters.h>
#include <chargefw/parameters/models/parameter_key.h>
#include <chargefw/parameters/models/parameter_set.h>
#include <chargefw/parameters/models/parameter_set_metadata.h>
#include <chargefw/parameters/models/parameter_view.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <snitch/snitch.hpp>
#include <span>
#include <stdexcept>
#include <vector>

namespace parameters = chargefw::parameters;
namespace features = chargefw::features;
namespace methods = chargefw::methods;
namespace core = chargefw::core;
namespace charges = chargefw::charges;

namespace {

auto make_parameter_set(const double hydrogen_width = 1.0, const double oxygen_width = 1.0)
    -> parameters::ParameterSet {
    return chargefw::test::make_sqe_ho_parameters("sqeqp", {.hydrogen_width = hydrogen_width,
                                                            .oxygen_width = oxygen_width,
                                                            .hydrogen_q0 = 0.25,
                                                            .oxygen_q0 = -0.5});
}

auto calculate_sqeqp(const core::Molecule& molecule, const parameters::ParameterSet& parameter_set,
                     const double target_charge,
                     const std::span<const methods::FixedPointSource> fixed_sources = {},
                     const std::vector<std::size_t>& atom_parameter_indices = {0, 1},
                     const std::vector<std::size_t>& bond_parameter_indices = {0})
    -> charges::AtomicCharges {
    const features::PreparedMolecule prepared{molecule};
    const features::ConformerFeatures geometry{molecule, 0};
    const parameters::ParameterClassification classification{
        parameters::AtomParameterClassification{atom_parameter_indices},
        parameters::BondParameterClassification{bond_parameter_indices}};
    const parameters::ParameterView parameter_view{parameter_set, classification};
    const auto* method = methods::method_registry().find("sqeqp");
    if (method == nullptr) {
        throw std::logic_error{"SQE+qp is not registered"};
    }
    const auto options = methods::make_default_options(method->option_schema());
    const methods::CalculationInput input{prepared,  options,         target_charge,
                                          &geometry, &parameter_view, fixed_sources};
    return method->calculate(input);
}

auto make_disconnected_parameters() -> parameters::ParameterSet {
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{
            .id = "disconnected-sqeqp", .method_id = "sqeqp", .name = "Disconnected SQE+qp"},
        {},
        parameters::AtomParameters{{{.key = chargefw::test::plain_atom_key(1),
                                     .parameters = {{.name = "electronegativity", .value = 1.0},
                                                    {.name = "hardness", .value = 5.0},
                                                    {.name = "width", .value = 1.0},
                                                    {.name = "q0", .value = 0.2}}},
                                    {.key = chargefw::test::plain_atom_key(8),
                                     .parameters = {{.name = "electronegativity", .value = 2.0},
                                                    {.name = "hardness", .value = 6.0},
                                                    {.name = "width", .value = 1.0},
                                                    {.name = "q0", .value = 0.1}}},
                                    {.key = chargefw::test::plain_atom_key(6),
                                     .parameters = {{.name = "electronegativity", .value = 3.0},
                                                    {.name = "hardness", .value = 7.0},
                                                    {.name = "width", .value = 0.8},
                                                    {.name = "q0", .value = -0.2}}},
                                    {.key = chargefw::test::plain_atom_key(7),
                                     .parameters = {{.name = "electronegativity", .value = 4.0},
                                                    {.name = "hardness", .value = 8.0},
                                                    {.name = "width", .value = 0.9},
                                                    {.name = "q0", .value = 0.1}}}}},
        parameters::BondParameters{{{.key = chargefw::test::single_bond_key(1, 8),
                                     .parameters = {{.name = "kappa", .value = 1.2}}},
                                    {.key = chargefw::test::single_bond_key(6, 7),
                                     .parameters = {{.name = "kappa", .value = 0.9}}}}}};
}

} // namespace

TEST_CASE("SQE+qp subtracts the Gaussian-to-point source field before bond projection",
          "[methods][sqeqp]") {
    const auto molecule = core::Molecule{std::vector{core::Atom{1}, core::Atom{8}},
                                         {core::Bond{0, 1}},
                                         {core::Conformer{{{0.0, 0.0, 0.0}, {1.5, 0.0, 0.0}}}},
                                         "two-atom-reference"};
    constexpr auto target_charge = 0.7;
    constexpr auto seed_h = 0.25 - ((0.25 - 0.5) - target_charge) / 2.0;
    constexpr auto seed_o = -0.5 - ((0.25 - 0.5) - target_charge) / 2.0;
    const auto interaction = [](const double distance, const double width_a, const double width_b) {
        const auto width_sum = 2.0 * width_a * width_a + 2.0 * width_b * width_b;
        return width_sum == 0.0 ? 1.0 / distance
                                : std::erf(distance / std::sqrt(width_sum)) / distance;
    };
    auto positive_width_result = std::vector<double>{};
    const auto source_cases = std::array{methods::FixedPointSource{{0.0, 2.0, 0.0}, 0.4},
                                         methods::FixedPointSource{{0.0, 2.0, 0.0}, -0.4},
                                         methods::FixedPointSource{{2.5, 0.0, 0.0}, 0.4},
                                         methods::FixedPointSource{{0.0, 2.0, 0.0}, 0.0}};
    const auto no_source = calculate_sqeqp(molecule, make_parameter_set(1.2, -0.65), target_charge);

    for (const auto hydrogen_width : {1.2, -1.2, 0.0}) {
        CAPTURE(hydrogen_width);
        const auto parameter_set = make_parameter_set(hydrogen_width, -0.65);
        const auto active_interaction = interaction(1.5, hydrogen_width, -0.65);
        const auto denominator = 13.8904 + 13.364 - 2.0 * active_interaction + 1.0;
        // Source variations and width limits are independent; only the reference needs all widths.
        const auto cases = std::span{source_cases}.first(hydrogen_width == 1.2 ? 4 : 1);

        for (const auto& source : cases) {
            CAPTURE(source.charge, source.position.x);
            const auto input_sources = std::span<const methods::FixedPointSource>{&source, 1};
            const auto charges =
                calculate_sqeqp(molecule, parameter_set, target_charge, input_sources);
            const auto distance_h = core::distance(core::Position{0.0, 0.0, 0.0}, source.position);
            const auto distance_o = core::distance(core::Position{1.5, 0.0, 0.0}, source.position);
            const auto potential_h = source.charge * interaction(distance_h, hydrogen_width, 0.0);
            const auto potential_o = source.charge * interaction(distance_o, -0.65, 0.0);
            const auto rhs_h = -4.5280 - potential_h - active_interaction * seed_o;
            const auto rhs_o = -8.741 - potential_o - active_interaction * seed_h;
            const auto transfer = (rhs_h - rhs_o) / denominator;
            const auto expected_h = seed_h + transfer;
            const auto expected_o = seed_o - transfer;

            CHECK(std::abs(charges[0] - expected_h) < 1.0e-12);
            CHECK(std::abs(charges[1] - expected_o) < 1.0e-12);
            CHECK(std::abs(charges.total() - target_charge) < 1.0e-12);
            if (source.charge == 0.0) {
                CHECK(std::ranges::equal(charges.values(), no_source.values()));
            }
            const auto is_positive_reference_source =
                source.charge == 0.4 && source.position.x == 0.0 && source.position.y == 2.0;
            if (is_positive_reference_source && hydrogen_width == 1.2) {
                positive_width_result.assign(charges.values().begin(), charges.values().end());
            } else if (is_positive_reference_source && hydrogen_width == -1.2) {
                CHECK(std::ranges::equal(charges.values(), positive_width_result));
            }
            if (source.charge == -0.4 || source.position.x == 2.5) {
                CHECK(std::abs(charges[0] - positive_width_result[0]) > 1.0e-5);
            }
        }
    }
}

TEST_CASE("SQE+qp keeps each disconnected component at its formal charge", "[methods][sqeqp]") {
    // Two bonded pairs 3 Å apart, with formal charges +1 and -1, polarized by a fixed source.
    const auto molecule = core::Molecule{
        std::vector{core::Atom{1}, core::Atom{8, 1}, core::Atom{6, -1}, core::Atom{7}},
        {core::Bond{0, 1}, core::Bond{2, 3}},
        {core::Conformer{{{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {4.0, 0.0, 0.0}, {5.0, 0.0, 0.0}}}},
        "two-components"};
    const auto parameter_set = make_disconnected_parameters();
    const auto sources = std::array{methods::FixedPointSource{{0.0, 2.0, 0.0}, 0.3}};

    for (const auto& fixed_sources : {std::span<const methods::FixedPointSource>{},
                                      std::span<const methods::FixedPointSource>{sources}}) {
        CAPTURE(fixed_sources.size());
        const auto charges =
            calculate_sqeqp(molecule, parameter_set, 0.0, fixed_sources, {0, 1, 2, 3}, {0, 1});
        CHECK(std::abs((charges[0] + charges[1]) - 1.0) < 1.0e-12);
        CHECK(std::abs((charges[2] + charges[3]) + 1.0) < 1.0e-12);
    }

    // A total charge other than the formal charges is spread over all atoms.
    const auto shifted = calculate_sqeqp(molecule, parameter_set, 0.8, {}, {0, 1, 2, 3}, {0, 1});
    CHECK(std::abs((shifted[0] + shifted[1]) - 1.4) < 1.0e-12);
    CHECK(std::abs((shifted[2] + shifted[3]) + 0.6) < 1.0e-12);
}

TEST_CASE("SQE+qp starts from given reference charges as they are", "[methods][sqeqp]") {
    const auto molecule = core::Molecule{std::vector{core::Atom{1}, core::Atom{8}},
                                         {core::Bond{0, 1}},
                                         {core::Conformer{{{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}}}},
                                         "bonded-pair"};
    const auto parameter_set = make_parameter_set();
    const features::PreparedMolecule prepared{molecule};
    const features::ConformerFeatures geometry{molecule, 0};
    const parameters::ParameterClassification classification{
        parameters::AtomParameterClassification{{0, 1}},
        parameters::BondParameterClassification{{0}}};
    const parameters::ParameterView parameter_view{parameter_set, classification};
    const auto* method = methods::method_registry().find("sqeqp");
    REQUIRE(method != nullptr);
    const auto options = methods::make_default_options(method->option_schema());

    // The result keeps the sum of the reference charges, not the target charge.
    const auto references = std::array{0.4, -0.1};
    const methods::CalculationInput input{prepared,        options, 0.0,       &geometry,
                                          &parameter_view, {},      references};
    CHECK(std::abs(method->calculate(input).total() - 0.3) < 1.0e-12);

    const auto too_few = std::array{0.0};
    CHECK_THROWS_AS((methods::CalculationInput{
                        prepared, options, 0.0, &geometry, &parameter_view, {}, too_few}),
                    std::invalid_argument);
}

TEST_CASE("SQE+qp validates sources before the no-bond return", "[methods][sqeqp]") {
    const auto molecule = core::Molecule{std::vector{core::Atom{1}, core::Atom{8}},
                                         {},
                                         {core::Conformer{{{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}}}},
                                         "isolated-pair"};
    const auto parameter_set = make_parameter_set();
    const auto invalid_charge = std::array{
        methods::FixedPointSource{{0.0, 2.0, 0.0}, std::numeric_limits<double>::quiet_NaN()}};
    CHECK_THROWS_AS(calculate_sqeqp(molecule, parameter_set, 0.7, invalid_charge, {0, 1}, {}),
                    std::invalid_argument);
    const auto invalid_position = std::array{
        methods::FixedPointSource{{std::numeric_limits<double>::quiet_NaN(), 2.0, 0.0}, 0.4}};
    CHECK_THROWS_AS(calculate_sqeqp(molecule, parameter_set, 0.7, invalid_position, {0, 1}, {}),
                    std::invalid_argument);

    const auto nonfinite_geometry = core::Molecule{
        std::vector{core::Atom{1}, core::Atom{8}},
        {},
        {core::Conformer{{{std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0}, {1.0, 0.0, 0.0}}}},
        "nonfinite-geometry"};
    const auto valid_source = std::array{methods::FixedPointSource{{0.0, 2.0, 0.0}, 0.4}};
    CHECK_THROWS_AS(
        calculate_sqeqp(nonfinite_geometry, parameter_set, 0.7, valid_source, {0, 1}, {}),
        std::invalid_argument);

    const auto empty_sources = std::span<const methods::FixedPointSource>{};
    const auto no_field = calculate_sqeqp(molecule, parameter_set, 0.7, empty_sources, {0, 1}, {});
    const auto with_source =
        calculate_sqeqp(molecule, parameter_set, 0.7, valid_source, {0, 1}, {});
    // Unbonded atoms are separate components at their formal charges, 0; the rest of the total
    // charge is spread evenly.
    CHECK(std::abs(no_field[0] - 0.35) < 1.0e-12);
    CHECK(std::abs(no_field[1] - 0.35) < 1.0e-12);
    CHECK(std::ranges::equal(no_field.values(), with_source.values()));
}
