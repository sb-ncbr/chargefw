#include "support/test_methods.h"
#include "support/test_molecules.h"

#include <chargefw/core/atom.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/core/position.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/calculation_input.h>
#include <chargefw/methods/method.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/methods/method_prerequisites.h>
#include <chargefw/methods/method_registry.h>

#include <limits>
#include <optional>
#include <snitch/snitch.hpp>
#include <span>
#include <string>
#include <utility>
#include <vector>

namespace methods = chargefw::methods;

TEST_CASE("prerequisite issue kinds convert to stable strings", "[methods][prerequisites]") {
    CHECK(methods::to_string(methods::PrerequisiteIssueKind::invalid_options) == "invalid_options");
    CHECK(methods::to_string(methods::PrerequisiteIssueKind::missing_parameters) ==
          "missing_parameters");
    CHECK(methods::to_string(methods::PrerequisiteIssueKind::unsupported_fixed_ions) ==
          "unsupported_fixed_ions");
}

namespace core = chargefw::core;
namespace features = chargefw::features;

namespace {

auto make_collection() -> core::MoleculeCollection {
    std::vector molecules{chargefw::test::make_water(),
                          chargefw::test::make_formally_charged_pair()};

    return core::MoleculeCollection{std::move(molecules), "test-collection"};
}

} // namespace

TEST_CASE("parameterless built-ins accept water", "[methods][method-prerequisites]") {
    const auto water = chargefw::test::make_water();
    const features::PreparedMolecule prepared_water{water};
    const auto empty_options = methods::MethodOptions{};

    for (const auto method_id : {"dummy", "formal"}) {
        CAPTURE(method_id);
        const auto* method = methods::method_registry().find(method_id);
        REQUIRE(method != nullptr);
        CHECK(method->check_method_prerequisites(
            {.prepared_molecule = prepared_water, .method_options = empty_options}));
    }
}

TEST_CASE("coordinate requirements detect molecules without conformers",
          "[methods][method-prerequisites]") {
    const chargefw::test::StubMethod coordinates_method{"coordinates-test", {.coordinates = true}};
    const auto empty_options = methods::MethodOptions{};
    const auto water = chargefw::test::make_water();
    const auto charged_pair = chargefw::test::make_formally_charged_pair();
    const features::PreparedMolecule prepared_water{water};
    const features::PreparedMolecule prepared_charged_pair{charged_pair};

    CHECK(coordinates_method.check_method_prerequisites(
        {.prepared_molecule = prepared_water, .method_options = empty_options}));

    const auto missing_coordinates = coordinates_method.check_method_prerequisites(
        {.prepared_molecule = prepared_charged_pair, .method_options = empty_options});
    CHECK(!missing_coordinates);
    REQUIRE(missing_coordinates.issues().size() == 1);
    CHECK(missing_coordinates.issues()[0].kind == methods::PrerequisiteIssueKind::missing_feature);

    const auto collection = make_collection();
    const features::PreparedMoleculeCollection prepared_collection{collection};
    const auto* dummy = methods::method_registry().find("dummy");
    REQUIRE(dummy != nullptr);
    CHECK(methods::check_method_prerequisites(*dummy, prepared_collection, empty_options));

    const auto collection_result =
        methods::check_method_prerequisites(coordinates_method, prepared_collection, empty_options);
    CHECK(!collection_result);
    REQUIRE(collection_result.issues().size() == 1);
    CHECK(collection_result.issues()[0].kind == methods::PrerequisiteIssueKind::missing_feature);
    CHECK(collection_result.issues()[0].molecule_index == 1);
    CHECK(collection_result.issues()[0].message.contains("molecule 2"));
}

TEST_CASE("bond requirements reject molecules with several atoms and no bonds",
          "[methods][method-prerequisites]") {
    const chargefw::test::StubMethod bonds_method{"bonds-test", {.bonds = true}};
    const auto empty_options = methods::MethodOptions{};
    const auto water = chargefw::test::make_water();
    const core::Molecule unbonded_pair{std::vector{core::Atom{11}, core::Atom{17}}};
    const core::Molecule single_atom{std::vector{core::Atom{11}}};
    const features::PreparedMolecule prepared_water{water};
    const features::PreparedMolecule prepared_unbonded_pair{unbonded_pair};
    const features::PreparedMolecule prepared_single_atom{single_atom};

    CHECK(bonds_method.check_method_prerequisites(
        {.prepared_molecule = prepared_water, .method_options = empty_options}));
    CHECK(bonds_method.check_method_prerequisites(
        {.prepared_molecule = prepared_single_atom, .method_options = empty_options}));

    const auto missing_bonds = bonds_method.check_method_prerequisites(
        {.prepared_molecule = prepared_unbonded_pair, .method_options = empty_options});
    CHECK(!missing_bonds);
    REQUIRE(missing_bonds.issues().size() == 1);
    CHECK(missing_bonds.issues()[0].kind == methods::PrerequisiteIssueKind::missing_feature);

    // Methods built on the bond graph declare the requirement; geometry-only methods do not.
    for (const auto* method_id : {"peoe", "delre", "sqeqp", "abeem"}) {
        CAPTURE(method_id);
        const auto* method = methods::method_registry().find(method_id);
        REQUIRE(method != nullptr);
        CHECK(method->requirements().bonds);
    }
    for (const auto* method_id : {"eem", "qeq", "veem"}) {
        CAPTURE(method_id);
        const auto* method = methods::method_registry().find(method_id);
        REQUIRE(method != nullptr);
        CHECK_FALSE(method->requirements().bonds);
    }
}

TEST_CASE("coordinate requirements reject coincident atoms in every conformer",
          "[methods][method-prerequisites]") {
    const chargefw::test::StubMethod coordinates_method{"coordinates-test", {.coordinates = true}};
    const core::Molecule coincident_atoms{
        std::vector{core::Atom{1}, core::Atom{1}},
        {},
        std::vector{core::Conformer{{core::Position{.x = 1.5, .y = -2.25, .z = 3.75},
                                     core::Position{.x = 1.5, .y = -2.25, .z = 3.75}},
                                    "first"},
                    core::Conformer{{core::Position{}, core::Position{.x = -0.0}}, "second"}}};
    const features::PreparedMolecule prepared{coincident_atoms};

    const auto result = coordinates_method.check_method_prerequisites(
        {.prepared_molecule = prepared, .method_options = methods::MethodOptions{}});
    CHECK(!result);
    REQUIRE(result.issues().size() == 2);
    for (std::size_t conformer_index = 0; conformer_index < 2; ++conformer_index) {
        const auto& issue = result.issues()[conformer_index];
        CHECK(issue.kind == methods::PrerequisiteIssueKind::invalid_geometry);
        CHECK(issue.atom_index == 1);
        CHECK(issue.conformer_index == conformer_index);
        CHECK(issue.message.contains("conformer " + std::to_string(conformer_index + 1)));
    }
    CHECK(result.issues()[0].message.contains("atom 1 (H, formal charge 0)"));
}

TEST_CASE("coordinate requirements reject non-finite coordinates",
          "[methods][method-prerequisites]") {
    const chargefw::test::StubMethod coordinates_method{"coordinates-test", {.coordinates = true}};
    const core::Molecule nonfinite_atom{
        std::vector{core::Atom{1}},
        {},
        std::vector{core::Conformer{{core::Position{.x = std::numeric_limits<double>::quiet_NaN()}},
                                    "nonfinite"}}};
    const features::PreparedMolecule prepared{nonfinite_atom};

    const auto result = coordinates_method.check_method_prerequisites(
        {.prepared_molecule = prepared, .method_options = methods::MethodOptions{}});
    CHECK(!result);
    REQUIRE(result.issues().size() == 1);
    const auto& issue = result.issues()[0];
    CHECK(issue.kind == methods::PrerequisiteIssueKind::invalid_geometry);
    CHECK(issue.atom_index == 0);
    CHECK(issue.conformer_index == 0);
    CHECK(issue.message.contains("method 'coordinates-test'"));
    CHECK(issue.message.contains("atom 1 (H, formal charge 0)"));
}

TEST_CASE("resource requirements do not restrict prerequisites",
          "[methods][method-prerequisites]") {
    const chargefw::test::StubMethod dense_method{
        "dense-test",
        {.resources = {.time = methods::ComplexityTerm::atoms_cubed,
                       .memory = methods::ComplexityTerm::atoms_squared}}};
    const auto water = chargefw::test::make_water();
    const features::PreparedMolecule prepared_water{water};

    const auto result = dense_method.check_method_prerequisites(
        {.prepared_molecule = prepared_water, .method_options = methods::MethodOptions{}});
    CHECK(result);
    CHECK(result.issues().empty());
}

TEST_CASE("element prerequisites identify unsupported atoms", "[methods][method-prerequisites]") {
    const auto* abeem = methods::method_registry().find("abeem");
    REQUIRE(abeem != nullptr);
    const core::Molecule berkelium_molecule{
        std::vector{core::Atom{97}}, {}, std::vector{core::Conformer{{core::Position{}}}}};
    const features::PreparedMolecule prepared_berkelium{berkelium_molecule};

    const auto result = abeem->check_method_prerequisites(
        {.prepared_molecule = prepared_berkelium, .method_options = methods::MethodOptions{}});
    CHECK(!result);
    REQUIRE(result.issues().size() == 1);
    CHECK(result.issues()[0].kind == methods::PrerequisiteIssueKind::unsupported_molecule);
    CHECK(result.issues()[0].atom_index == 0);
    CHECK(result.issues()[0].message.contains("atom 1 (Bk, formal charge 0)"));
}
