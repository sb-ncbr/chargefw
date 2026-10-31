#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/molecule.h>
#include <chargefw/features/topology_features.h>
#include <chargefw/parameters/classification/parameter_classifier.h>
#include <chargefw/parameters/models/atom_parameters.h>
#include <chargefw/parameters/models/bond_parameters.h>
#include <chargefw/parameters/models/parameter_key.h>
#include <chargefw/parameters/models/parameter_set.h>
#include <chargefw/parameters/models/parameter_set_metadata.h>

#include <string>
#include <utility>
#include <vector>

#include <snitch/snitch.hpp>

namespace core = chargefw::core;
namespace features = chargefw::features;
namespace parameters = chargefw::parameters;

namespace {

auto test_metadata() -> parameters::ParameterSetMetadata {
    return {.id = "test-parameters", .method_id = "test-method", .name = "Test parameters"};
}

auto bond_type_key(const parameters::BondParameterClassificationKind classification,
                   std::string type) -> parameters::BondTypeKey {
    return {.classification = classification, .type = std::move(type)};
}

auto classify(const core::Molecule& molecule, const parameters::ParameterSet& parameter_set)
    -> parameters::ParameterClassification {
    const features::TopologyFeatures topology{molecule};

    return parameters::classify_parameters(molecule, topology, parameter_set);
}

auto make_c_cl_f_h_fragment() -> core::Molecule {
    std::vector atoms{core::Atom{6, 0, "C"}, core::Atom{17, 0, "Cl"}, core::Atom{9, 0, "F"},
                      core::Atom{1, 0, "H"}};

    std::vector bonds{core::Bond{0, 1, core::BondOrder::SINGLE},
                      core::Bond{0, 2, core::BondOrder::SINGLE},
                      core::Bond{0, 3, core::BondOrder::SINGLE}};

    return core::Molecule{std::move(atoms), std::move(bonds), {}, "c-cl-f-h"};
}

auto make_carbonyl() -> core::Molecule {
    return core::Molecule{{core::Atom{6, 0, "C"}, core::Atom{8, 0, "O"}},
                          {core::Bond{0, 1, core::BondOrder::DOUBLE}},
                          {},
                          "carbonyl"};
}

auto make_lone_oxygen() -> core::Molecule {
    return core::Molecule{{core::Atom{8, 0, "O"}}, {}, {}, "lone-oxygen"};
}

} // namespace

TEST_CASE("parameter classifier maps water atoms and bonds to parameter entries",
          "[parameters][classifier]") {
    const auto water = chargefw::test::make_water();

    const parameters::ParameterSet water_parameters{
        test_metadata(),
        {},
        parameters::AtomParameters{
            {{.key = chargefw::test::atom_key(
                  1, parameters::AtomParameterClassificationKind::BONDED_ELEMENTS, "O"),
              .parameters = {{.name = "value", .value = 1.0}}},
             {.key = chargefw::test::atom_key(
                  8, parameters::AtomParameterClassificationKind::BONDED_ELEMENTS, "HH"),
              .parameters = {{.name = "value", .value = 2.0}}}}},
        parameters::BondParameters{
            {{.key = {.first_atom = chargefw::test::atom_key(
                          1, parameters::AtomParameterClassificationKind::BONDED_ELEMENTS, "O"),
                      .second_atom = chargefw::test::atom_key(
                          8, parameters::AtomParameterClassificationKind::BONDED_ELEMENTS, "HH"),
                      .bond = bond_type_key(parameters::BondParameterClassificationKind::BOND_ORDER,
                                            "1")},
              .parameters = {{.name = "value", .value = 3.0}}}}}};

    const auto water_classification = classify(water, water_parameters);

    REQUIRE(water_classification.atom().size() == 3);
    CHECK(water_classification.atom()[0] == 1);
    CHECK(water_classification.atom()[1] == 0);
    CHECK(water_classification.atom()[2] == 0);

    REQUIRE(water_classification.bond().size() == 2);
    CHECK(water_classification.bond()[0] == 0);
    CHECK(water_classification.bond()[1] == 0);
}

TEST_CASE("parameter classifier maps mixed element fragment atoms", "[parameters][classifier]") {
    const auto mixed = make_c_cl_f_h_fragment();

    const parameters::ParameterSet mixed_parameters{
        test_metadata(),
        {},
        parameters::AtomParameters{
            {{.key = chargefw::test::atom_key(
                  6, parameters::AtomParameterClassificationKind::BONDED_ELEMENTS, "ClFH"),
              .parameters = {{.name = "value", .value = 1.0}}},
             {.key = chargefw::test::atom_key(
                  17, parameters::AtomParameterClassificationKind::BONDED_ELEMENTS, "C"),
              .parameters = {{.name = "value", .value = 2.0}}},
             {.key = chargefw::test::atom_key(
                  9, parameters::AtomParameterClassificationKind::BONDED_ELEMENTS, "C"),
              .parameters = {{.name = "value", .value = 3.0}}},
             {.key = chargefw::test::atom_key(
                  1, parameters::AtomParameterClassificationKind::BONDED_ELEMENTS, "C"),
              .parameters = {{.name = "value", .value = 4.0}}}}}};

    const auto mixed_classification = classify(mixed, mixed_parameters);

    REQUIRE(mixed_classification.atom().size() == 4);
    CHECK(mixed_classification.atom()[0] == 0);
    CHECK(mixed_classification.atom()[1] == 1);
    CHECK(mixed_classification.atom()[2] == 2);
    CHECK(mixed_classification.atom()[3] == 3);
}

TEST_CASE("parameter classifier uses permissive bond-order fallback only when requested",
          "[parameters][classifier]") {
    const auto carbonyl = make_carbonyl();
    const parameters::ParameterSet parameter_set{
        test_metadata(),
        {},
        parameters::AtomParameters{{{.key = chargefw::test::hbo_atom_key(6, "1"),
                                     .parameters = {{.name = "value", .value = 1.0}}},
                                    {.key = chargefw::test::hbo_atom_key(8, "1"),
                                     .parameters = {{.name = "value", .value = 2.0}}}}},
        parameters::BondParameters{{{.key = chargefw::test::single_bond_key(6, 8),
                                     .parameters = {{.name = "value", .value = 3.0}}}}}};

    const auto exact = parameters::try_classify_parameters(
        carbonyl, features::TopologyFeatures{carbonyl}, parameter_set);
    CHECK_FALSE(exact);
    REQUIRE(exact.issues().size() == 3);
    CHECK(exact.issues()[0].kind == parameters::ClassificationIssueKind::MISSING_ATOM_PARAMETER);
    CHECK(exact.issues()[2].kind == parameters::ClassificationIssueKind::MISSING_BOND_PARAMETER);

    const auto permissive = parameters::try_classify_parameters(
        carbonyl, features::TopologyFeatures{carbonyl}, parameter_set, {.permissive_types = true});
    REQUIRE(permissive);
    const auto atom_indices = permissive.classification().atom().parameter_entry_indices();
    const auto bond_indices = permissive.classification().bond().parameter_entry_indices();
    REQUIRE(atom_indices.size() == 2);
    CHECK(atom_indices[0] == 0);
    CHECK(atom_indices[1] == 1);
    REQUIRE(bond_indices.size() == 1);
    CHECK(bond_indices[0] == 0);

    // An atom without bonds has highest bond order 0 and falls back to the single-bond entry.
    const auto lone_oxygen = make_lone_oxygen();
    const features::TopologyFeatures lone_oxygen_topology{lone_oxygen};
    CHECK_FALSE(
        parameters::try_classify_parameters(lone_oxygen, lone_oxygen_topology, parameter_set));
    const auto lone_permissive = parameters::try_classify_parameters(
        lone_oxygen, lone_oxygen_topology, parameter_set, {.permissive_types = true});
    REQUIRE(lone_permissive);
    CHECK(lone_permissive.classification().atom()[0] == 1);
}

TEST_CASE("parameter classifier resolves wildcard entries in declaration order",
          "[parameters][classifier]") {
    const auto water = chargefw::test::make_water_graph();
    const parameters::ParameterSet parameter_set{
        test_metadata(),
        {},
        parameters::AtomParameters{{{.key = chargefw::test::plain_atom_key(0),
                                     .parameters = {{.name = "value", .value = 1.0}}},
                                    {.key = chargefw::test::plain_atom_key(8),
                                     .parameters = {{.name = "value", .value = 2.0}}}}}};

    const auto classification = classify(water, parameter_set);
    const auto atom_indices = classification.atom().parameter_entry_indices();
    REQUIRE(atom_indices.size() == 3);
    CHECK(atom_indices[0] == 0);
    CHECK(atom_indices[1] == 0);
    CHECK(atom_indices[2] == 0);
}

TEST_CASE("parameter classifier reports missing entries with source indices",
          "[parameters][classifier]") {
    const auto water = chargefw::test::make_water_graph();
    const parameters::ParameterSet parameter_set{
        test_metadata(),
        {},
        parameters::AtomParameters{{{.key = chargefw::test::plain_atom_key(1),
                                     .parameters = {{.name = "value", .value = 1.0}}}}}};

    const auto result = parameters::try_classify_parameters(
        water, features::TopologyFeatures{water}, parameter_set);
    CHECK_FALSE(result);
    REQUIRE(result.issues().size() == 1);
    CHECK(result.issues()[0].kind == parameters::ClassificationIssueKind::MISSING_ATOM_PARAMETER);
    CHECK(result.issues()[0].object_index == 0);
}
