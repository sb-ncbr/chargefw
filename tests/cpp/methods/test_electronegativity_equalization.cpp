#include "support/test_calculation.h"
#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <chargefw/core/position.h>
#include <chargefw/features/conformer_features.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/methods/calculation_input.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/methods/method_registry.h>
#include <chargefw/parameters/classification/parameter_classification.h>
#include <chargefw/parameters/models/atom_parameters.h>
#include <chargefw/parameters/models/common_parameters.h>
#include <chargefw/parameters/models/parameter_key.h>
#include <chargefw/parameters/models/parameter_set.h>
#include <chargefw/parameters/models/parameter_set_metadata.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <optional>
#include <snitch/snitch.hpp>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace parameters = chargefw::parameters;
namespace core = chargefw::core;
namespace features = chargefw::features;
namespace methods = chargefw::methods;

namespace {

auto make_active_hydrogen_oxygen(const core::Position hydrogen = {},
                                 const core::Position oxygen = core::Position{.x = 2.0})
    -> core::Molecule {
    return core::Molecule{
        {core::Atom{1, 1}, core::Atom{8}}, {}, {core::Conformer{{hydrogen, oxygen}}}};
}

auto calculate_eem(const core::Molecule& molecule,
                   const std::span<const methods::FixedPointSource> fixed_sources,
                   const double target_charge, const double kappa = 1.7)
    -> chargefw::charges::AtomicCharges {
    const auto parameter_set = chargefw::test::make_eem_ho_parameters(kappa);
    const auto classification =
        parameters::ParameterClassification{parameters::AtomParameterClassification{
            molecule.atom_count() == 0 ? std::vector<std::size_t>{}
                                       : std::vector<std::size_t>{0, 1}}};
    const auto parameter_view = parameters::ParameterView{parameter_set, classification};
    const auto prepared = features::PreparedMolecule{molecule};
    auto geometry = std::optional<features::ConformerFeatures>{};
    if (molecule.conformer_count() != 0) {
        geometry.emplace(molecule);
    }
    const auto options = methods::MethodOptions{};
    const auto input = methods::CalculationInput{
        prepared,        options,      target_charge, geometry.has_value() ? &*geometry : nullptr,
        &parameter_view, fixed_sources};
    const auto* method = methods::method_registry().find("eem");
    if (method == nullptr) {
        throw std::runtime_error{"EEM method is missing from the builtin registry"};
    }
    return method->calculate(input);
}

} // namespace

TEST_CASE("EEM enforces a charged molecular target", "[methods][eem]") {
    const chargefw::core::Molecule cation{
        {chargefw::core::Atom{1, 1}},
        {},
        {chargefw::core::Conformer{{chargefw::core::Position{}}}}};
    const auto charge_set = chargefw::test::calculate_single_method(
        cation, "eem", {chargefw::test::make_eem_ho_parameters(1.0)});

    CHECK(std::abs(charge_set.assignment(0).charges[0] - 1.0) < 1.0e-12);
    CHECK(std::abs(charge_set.assignment(0).charges.total() - 1.0) < 1.0e-12);
}

TEST_CASE("EEM two-atom charges obey equalization at different distances and targets",
          "[methods][eem]") {
    for (const int target : {0, 1}) {
        CAPTURE(target);
        const chargefw::core::Molecule molecule{
            {chargefw::core::Atom{1, target}, chargefw::core::Atom{8}},
            {},
            {chargefw::core::Conformer{
                 {chargefw::core::Position{}, chargefw::core::Position{.x = 1.0}}},
             chargefw::core::Conformer{
                 {chargefw::core::Position{}, chargefw::core::Position{.x = 2.0}}}}};
        const auto charge_set = chargefw::test::calculate_method(
            molecule, "eem", {chargefw::test::make_eem_ho_parameters(1.0)});
        chargefw::test::assert_conformer_dependent(charge_set, 2);

        // EEM uses A = chi, B = hardness (no factor of two), and J = kappa/r.
        // The equalization equations give q_H = (1 + (9-J)*Q)/(14-2*J) and q_O = Q-q_H.
        // At r=1: Q=0 -> (1/12,-1/12), Q=1 -> (3/4,1/4).
        // At r=2: Q=0 -> (1/13,-1/13), Q=1 -> (19/26,7/26).
        const std::vector expected_hydrogen{target == 0 ? 1.0 / 12.0 : 3.0 / 4.0,
                                            target == 0 ? 1.0 / 13.0 : 19.0 / 26.0};
        for (std::size_t index = 0; index < expected_hydrogen.size(); ++index) {
            CAPTURE(index);
            const auto& charges = charge_set.assignment(index).charges;
            REQUIRE(charges.size() == 2);
            CHECK(std::abs(charges[0] - expected_hydrogen[index]) < 1.0e-12);
            CHECK(std::abs(charges[1] - (target - expected_hydrogen[index])) < 1.0e-12);
            CHECK(std::abs(charges.total() - target) < 1.0e-12);
        }
    }
}

TEST_CASE("EEM includes fixed point potentials in the constrained active solve", "[methods][eem]") {
    constexpr auto kappa = 1.7;
    constexpr auto target_charge = -0.35;
    const auto molecule = make_active_hydrogen_oxygen();
    const auto sources = std::array{methods::FixedPointSource{{-1.0, 0.0, 0.0}, 0.4},
                                    methods::FixedPointSource{{4.0, 1.0, 0.0}, -0.2}};
    const auto empty_span = std::span<const methods::FixedPointSource>{};
    const auto baseline = calculate_eem(molecule, empty_span, target_charge, kappa);
    const auto fixed = calculate_eem(molecule, std::span<const methods::FixedPointSource>{sources},
                                     target_charge, kappa);

    const auto J = kappa / 2.0;
    const auto phi_hydrogen =
        kappa * sources[0].charge / 1.0 + kappa * sources[1].charge / std::hypot(4.0, 1.0, 0.0);
    const auto phi_oxygen =
        kappa * sources[0].charge / 3.0 + kappa * sources[1].charge / std::hypot(2.0, 1.0, 0.0);
    const auto denominator = 5.0 + 9.0 - 2.0 * J;
    const auto expected_baseline_hydrogen = (2.0 - 1.0 + (9.0 - J) * target_charge) / denominator;
    const auto expected_fixed_hydrogen =
        (2.0 - 1.0 + phi_oxygen - phi_hydrogen + (9.0 - J) * target_charge) / denominator;

    REQUIRE(fixed.size() == 2);
    CHECK(std::abs(baseline[0] - expected_baseline_hydrogen) < 1.0e-12);
    CHECK(std::abs(fixed[0] - expected_fixed_hydrogen) < 1.0e-12);
    CHECK(std::abs(fixed[1] - (target_charge - expected_fixed_hydrogen)) < 1.0e-12);
    CHECK(fixed[0] < baseline[0]);
    CHECK(std::abs(fixed.total() - target_charge) < 1.0e-12);

    const auto zero_source = std::array{methods::FixedPointSource{{10.0, 2.0, 0.0}, 0.0}};
    const auto zero_source_charges = calculate_eem(
        molecule, std::span<const methods::FixedPointSource>{zero_source}, target_charge, kappa);
    CHECK(std::abs(zero_source_charges[0] - baseline[0]) < 1.0e-14);
    CHECK(std::abs(zero_source_charges[1] - baseline[1]) < 1.0e-14);
}

TEST_CASE("EEM rejects invalid fixed point fields and budgets", "[methods][eem]") {
    const auto molecule = make_active_hydrogen_oxygen();
    const auto check_invalid = [&](const core::Molecule& active_molecule,
                                   const std::span<const methods::FixedPointSource> sources,
                                   const double target_charge, const std::string_view diagnostic) {
        CHECK_THROWS_MATCHES(calculate_eem(active_molecule, sources, target_charge),
                             std::invalid_argument,
                             snitch::matchers::with_what_contains{diagnostic});
    };
    const auto nan = std::numeric_limits<double>::quiet_NaN();
    const auto infinity = std::numeric_limits<double>::infinity();
    const auto invalid_sources = std::array{methods::FixedPointSource{{4.0, 1.0, 0.0}, nan},
                                            methods::FixedPointSource{{nan, 1.0, 0.0}, 0.25},
                                            methods::FixedPointSource{{4.0, infinity, 0.0}, 0.25},
                                            methods::FixedPointSource{{4.0, 1.0, -infinity}, 0.25}};
    for (const auto& source : invalid_sources) {
        check_invalid(molecule, std::span{&source, 1}, -0.35, "fixed source 0");
    }

    const auto valid_source = std::array{methods::FixedPointSource{{4.0, 1.0, 0.0}, 0.25}};
    check_invalid(molecule, valid_source, std::numeric_limits<double>::infinity(),
                  "finite target charge");

    const auto empty_molecule = core::Molecule{std::vector<core::Atom>{}};
    check_invalid(empty_molecule, valid_source, -0.35, "requires active atoms");

    const auto nonfinite_active =
        make_active_hydrogen_oxygen({}, {std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0});
    check_invalid(nonfinite_active, valid_source, -0.35, "active atom 1");
}

TEST_CASE("EEM fixed-source validation requires geometry for its active molecule",
          "[methods][eem]") {
    const auto molecule = make_active_hydrogen_oxygen();
    const auto other_molecule = make_active_hydrogen_oxygen({}, {3.0, 0.0, 0.0});
    const auto parameter_set = chargefw::test::make_eem_ho_parameters(1.7);
    const auto classification = parameters::ParameterClassification{
        parameters::AtomParameterClassification{std::vector<std::size_t>{0, 1}}};
    const auto parameter_view = parameters::ParameterView{parameter_set, classification};
    const auto prepared = features::PreparedMolecule{molecule};
    const auto foreign_geometry = features::ConformerFeatures{other_molecule};
    const auto options = methods::MethodOptions{};
    const auto sources = std::array{methods::FixedPointSource{{4.0, 1.0, 0.0}, 0.25}};
    const auto source_span = std::span<const methods::FixedPointSource>{sources};
    const auto* eem = methods::method_registry().find("eem");
    REQUIRE(eem != nullptr);

    const auto missing_geometry =
        methods::CalculationInput{prepared, options, -0.35, nullptr, &parameter_view, source_span};
    CHECK_THROWS_AS(eem->calculate(missing_geometry), std::logic_error);

    const auto wrong_geometry = methods::CalculationInput{
        prepared, options, -0.35, &foreign_geometry, &parameter_view, source_span};
    CHECK_THROWS_MATCHES(eem->calculate(wrong_geometry), std::invalid_argument,
                         snitch::matchers::with_what_contains{
                             "fixed-source geometry does not belong to the active molecule"});
}

TEST_CASE("QEq overlap term changes charges while preserving neutrality", "[methods][qeq]") {
    const auto water = chargefw::test::make_water();
    const auto parameters = chargefw::test::make_qeq_ho_parameters();
    auto options = methods::MethodOptions{};
    options.set("overlap_term", std::string{"Nishimoto-Mataga"});
    const auto nishimoto_mataga =
        chargefw::test::calculate_single_method(water, "qeq", {parameters}, &options);
    options.set("overlap_term", std::string{"Ohno"});
    const auto ohno = chargefw::test::calculate_single_method(water, "qeq", {parameters}, &options);

    const auto& first = nishimoto_mataga.assignment(0).charges;
    const auto& second = ohno.assignment(0).charges;
    chargefw::test::assert_neutral_water_charges(first, 1.0e-10);
    chargefw::test::assert_neutral_water_charges(second, 1.0e-10);
    CHECK(std::abs(first[0] - second[0]) > 1.0e-6);
}
