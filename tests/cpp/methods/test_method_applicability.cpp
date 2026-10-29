#include "support/test_methods.h"
#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <chargefw/charges/atomic_charges.h>
#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method.h>
#include <chargefw/methods/method_applicability.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/methods/method_registry.h>
#include <chargefw/parameters/models/atom_parameters.h>
#include <chargefw/parameters/models/parameter_key.h>
#include <chargefw/parameters/models/parameter_set.h>
#include <chargefw/parameters/models/parameter_set_metadata.h>

#include <optional>
#include <snitch/snitch.hpp>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace core = chargefw::core;
namespace calculation = chargefw::calculation;
namespace features = chargefw::features;
namespace methods = chargefw::methods;
namespace parameters = chargefw::parameters;

TEST_CASE("execution applicability enums convert to stable strings", "[methods][applicability]") {
    CHECK(methods::to_string(methods::ExecutionAvailability::available_with_warning) ==
          "available_with_warning");
    CHECK(methods::to_string(methods::ExecutionIssueKind::unsupported_execution_mode) ==
          "unsupported_execution_mode");
}

namespace {

auto assessment_for(const methods::ApplicableMethod& candidate,
                    const calculation::ExecutionMode mode) -> const methods::ExecutionAssessment& {
    for (const auto& assessment : candidate.execution_assessments) {
        if (assessment.mode == mode) {
            return assessment;
        }
    }

    throw std::runtime_error{"missing execution assessment"};
}

auto make_collection() -> core::MoleculeCollection {
    std::vector molecules{chargefw::test::make_water(),
                          chargefw::test::make_formally_charged_pair()};

    return core::MoleculeCollection{std::move(molecules), "test-collection"};
}

auto make_wrong_method_parameters() -> parameters::ParameterSet {
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{.id = "wrong-method-parameters",
                                         .method_id = "other-method",
                                         .name = "Wrong method parameters"},
        {},
        parameters::AtomParameters{
            {{.key = chargefw::test::atom_key(1, parameters::AtomParameterClassificationKind::PLAIN,
                                              "*"),
              .parameters = {{.name = "value", .value = 1.0}}}}}};
}

auto make_water_only_parameters() -> parameters::ParameterSet {
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{.id = "water-only-parameters",
                                         .method_id = "atom-parameter-test",
                                         .name = "Water only parameters"},
        {},
        parameters::AtomParameters{
            {{.key = chargefw::test::atom_key(
                  1, parameters::AtomParameterClassificationKind::BONDED_ELEMENTS, "O"),
              .parameters = {{.name = "value", .value = 1.0}}},
             {.key = chargefw::test::atom_key(
                  8, parameters::AtomParameterClassificationKind::BONDED_ELEMENTS, "HH"),
              .parameters = {{.name = "value", .value = 2.0}}}}}};
}

auto make_collection_parameters() -> parameters::ParameterSet {
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{.id = "collection-parameters",
                                         .method_id = "atom-parameter-test",
                                         .name = "Collection parameters"},
        {},
        parameters::AtomParameters{
            {{.key = chargefw::test::atom_key(
                  1, parameters::AtomParameterClassificationKind::BONDED_ELEMENTS, "O"),
              .parameters = {{.name = "value", .value = 1.0}}},
             {.key = chargefw::test::atom_key(
                  8, parameters::AtomParameterClassificationKind::BONDED_ELEMENTS, "HH"),
              .parameters = {{.name = "value", .value = 2.0}}},
             {.key = chargefw::test::atom_key(7, parameters::AtomParameterClassificationKind::PLAIN,
                                              "*"),
              .parameters = {{.name = "value", .value = 3.0}}},
             {.key = chargefw::test::atom_key(
                  17, parameters::AtomParameterClassificationKind::PLAIN, "*"),
              .parameters = {{.name = "value", .value = 4.0}}}}}};
}

auto make_double_bonded_carbons() -> core::Molecule {
    return core::Molecule{std::vector{core::Atom{6}, core::Atom{6}},
                          std::vector{core::Bond{0, 1, core::BondOrder::DOUBLE}},
                          {},
                          "double-bonded-carbons"};
}

auto make_permissive_parameters(const bool include_exact_match) -> parameters::ParameterSet {
    auto entries = std::vector<parameters::AtomParameterEntry>{};

    if (include_exact_match) {
        entries.push_back(
            {.key = chargefw::test::atom_key(
                 6, parameters::AtomParameterClassificationKind::HIGHEST_BOND_ORDER, "2"),
             .parameters = {{.name = "value", .value = 2.0}}});
    }

    entries.push_back({.key = chargefw::test::atom_key(
                           6, parameters::AtomParameterClassificationKind::HIGHEST_BOND_ORDER, "1"),
                       .parameters = {{.name = "value", .value = 1.0}}});

    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{.id = "permissive-parameters",
                                         .method_id = "atom-parameter-test",
                                         .name = "Permissive parameters"},
        {},
        parameters::AtomParameters{std::move(entries)}};
}

} // namespace

TEST_CASE("applicability rejects invalid method options", "[methods][method-applicability]") {
    const auto collection = make_collection();
    const features::PreparedMoleculeCollection prepared_collection{collection};

    const auto& registry = methods::method_registry();
    const auto* gdac = registry.find("gdac");

    REQUIRE(gdac != nullptr);

    auto invalid_gdac_options = methods::MethodOptions{};
    invalid_gdac_options.set("iters", 0);
    const std::vector<const methods::Method*> gdac_methods{gdac};
    const auto invalid_gdac_result = methods::find_applicable_methods(
        {.molecules = prepared_collection,
         .methods = gdac_methods,
         .parameter_sets = {},
         .method_options = {{"gdac", std::move(invalid_gdac_options)}}});
    CHECK(invalid_gdac_result.applicable.empty());
    REQUIRE(invalid_gdac_result.rejected.size() == 1);
    REQUIRE(!invalid_gdac_result.rejected[0].issues.empty());
    CHECK(invalid_gdac_result.rejected[0].issues[0].kind ==
          methods::PrerequisiteIssueKind::invalid_options);
}

TEST_CASE("applicability pairs methods with compatible parameter sets",
          "[methods][method-applicability]") {
    const auto collection = make_collection();
    const features::PreparedMoleculeCollection prepared_collection{collection};
    const auto* dummy = methods::method_registry().find("dummy");
    REQUIRE(dummy != nullptr);

    const chargefw::test::StubMethod atom_parameter_method{"atom-parameter-test",
                                                           {.atom_parameters = {"value"}}};

    const std::vector<const methods::Method*> candidate_methods{dummy, &atom_parameter_method};

    const std::vector parameter_sets{make_wrong_method_parameters(), make_water_only_parameters(),
                                     make_collection_parameters()};

    const auto result = methods::find_applicable_methods({.molecules = prepared_collection,
                                                          .methods = candidate_methods,
                                                          .parameter_sets = parameter_sets});

    CHECK(!result.applicable.empty());
    REQUIRE(result.applicable.size() == 2);

    const auto& dummy_applicable = result.applicable[0];

    CHECK(dummy_applicable.method == dummy);
    CHECK(dummy_applicable.parameter_set == nullptr);
    CHECK(!dummy_applicable.uses_parameters());
    CHECK(dummy_applicable.classifications.empty());
    REQUIRE(dummy_applicable.execution_assessments.size() == 3);
    CHECK(assessment_for(dummy_applicable, calculation::ExecutionMode::full).availability ==
          methods::ExecutionAvailability::available);
    CHECK(assessment_for(dummy_applicable, calculation::ExecutionMode::cutoff).availability ==
          methods::ExecutionAvailability::unsupported);
    CHECK(assessment_for(dummy_applicable, calculation::ExecutionMode::cover).availability ==
          methods::ExecutionAvailability::unsupported);

    const auto& parameterized_applicable = result.applicable[1];

    CHECK(parameterized_applicable.method == &atom_parameter_method);
    REQUIRE(parameterized_applicable.parameter_set != nullptr);
    CHECK(parameterized_applicable.uses_parameters());
    CHECK(parameterized_applicable.parameter_set->id() ==
          std::string_view{"collection-parameters"});
    REQUIRE(parameterized_applicable.classifications.size() == prepared_collection.size());

    const auto& water_classification = parameterized_applicable.classifications[0];

    REQUIRE(water_classification.atom().size() == 3);
    CHECK(water_classification.atom()[0] == 1);
    CHECK(water_classification.atom()[1] == 0);
    CHECK(water_classification.atom()[2] == 0);

    const auto& charged_pair_classification = parameterized_applicable.classifications[1];

    REQUIRE(charged_pair_classification.atom().size() == 2);
    CHECK(charged_pair_classification.atom()[0] == 2);
    CHECK(charged_pair_classification.atom()[1] == 3);

    REQUIRE(result.rejected.size() == 1);

    CHECK(result.rejected[0].method_index == 1);
    CHECK(result.rejected[0].parameter_set_index == std::optional<std::size_t>{1});
    REQUIRE(!result.rejected[0].issues.empty());
    CHECK(result.rejected[0].issues[0].kind ==
          methods::PrerequisiteIssueKind::parameter_classification_failed);
    CHECK(result.rejected[0].issues[0].molecule_index == 1);
    CHECK(result.rejected[0].issues[0].atom_index == 0);
    CHECK(result.rejected[0].issues[0].message.contains("molecule 2"));
    CHECK(result.rejected[0].issues[0].message.contains(
        "atom 1 (source name 'N', N, formal charge +1)"));
    CHECK(result.rejected[0].issues[0].message.contains("water-only-parameters"));
}

TEST_CASE("applicability refuses methods that declare only common parameters",
          "[methods][method-applicability]") {
    const auto collection = make_collection();
    const features::PreparedMoleculeCollection prepared_collection{collection};
    const auto request = [&](const methods::Method& method) {
        const auto candidates = std::vector<const methods::Method*>{&method};
        return methods::find_applicable_methods(
            {.molecules = prepared_collection, .methods = candidates, .parameter_sets = {}});
    };

    // Common parameters fitted together with atom parameters belong in a parameter set.
    const chargefw::test::StubMethod with_atoms{
        "common-and-atom-test", {.common_parameters = {"k"}, .atom_parameters = {"value"}}};
    CHECK_NOTHROW(request(with_atoms));
    // Constants on their own are settings, so the method must declare them as options.
    const chargefw::test::StubMethod common_only{"common-only-test", {.common_parameters = {"k"}}};
    CHECK_THROWS_AS(request(common_only), std::invalid_argument);
}

TEST_CASE("applicability honors permissive parameter classification",
          "[methods][method-applicability]") {
    const chargefw::test::StubMethod atom_parameter_method{"atom-parameter-test",
                                                           {.atom_parameters = {"value"}}};

    const core::MoleculeCollection double_bonded_collection{
        std::vector{make_double_bonded_carbons()}};
    const features::PreparedMoleculeCollection double_bonded_prepared{double_bonded_collection};
    const std::vector<const methods::Method*> permissive_methods{&atom_parameter_method};
    const std::vector permissive_only_parameters{make_permissive_parameters(false)};

    const auto strict_result =
        methods::find_applicable_methods({.molecules = double_bonded_prepared,
                                          .methods = permissive_methods,
                                          .parameter_sets = permissive_only_parameters});
    CHECK(strict_result.applicable.empty());
    REQUIRE(strict_result.rejected.size() == 1);
    REQUIRE(!strict_result.rejected[0].issues.empty());
    CHECK(strict_result.rejected[0].issues[0].kind ==
          methods::PrerequisiteIssueKind::parameter_classification_failed);

    const auto permissive_result =
        methods::find_applicable_methods({.molecules = double_bonded_prepared,
                                          .methods = permissive_methods,
                                          .parameter_sets = permissive_only_parameters,
                                          .classification_options = {.permissive_types = true}});
    REQUIRE(permissive_result.applicable.size() == 1);
    REQUIRE(permissive_result.applicable[0].classifications.size() == 1);
    REQUIRE(permissive_result.applicable[0].classifications[0].atom().size() == 2);
    CHECK(permissive_result.applicable[0].classifications[0].atom()[0] == 0);
    CHECK(permissive_result.applicable[0].classifications[0].atom()[1] == 0);
    CHECK(assessment_for(permissive_result.applicable[0], calculation::ExecutionMode::full)
              .availability == methods::ExecutionAvailability::available);

    const std::vector exact_and_permissive_parameters{make_permissive_parameters(true)};
    const auto exact_result =
        methods::find_applicable_methods({.molecules = double_bonded_prepared,
                                          .methods = permissive_methods,
                                          .parameter_sets = exact_and_permissive_parameters,
                                          .classification_options = {.permissive_types = true}});
    REQUIRE(exact_result.applicable.size() == 1);
    REQUIRE(exact_result.applicable[0].classifications.size() == 1);
    REQUIRE(exact_result.applicable[0].classifications[0].atom().size() == 2);
    CHECK(exact_result.applicable[0].classifications[0].atom()[0] == 0);
    CHECK(exact_result.applicable[0].classifications[0].atom()[1] == 0);
    CHECK(
        assessment_for(exact_result.applicable[0], calculation::ExecutionMode::full).availability ==
        assessment_for(permissive_result.applicable[0], calculation::ExecutionMode::full)
            .availability);
}

TEST_CASE("applicability reports full-execution resource threshold warnings",
          "[methods][method-applicability]") {
    const auto collection = make_collection();
    const features::PreparedMoleculeCollection prepared_collection{collection};
    const auto* mgc = methods::method_registry().find("mgc");
    REQUIRE(mgc != nullptr);

    const chargefw::test::StubMethod inexpensive_method{"resource-test"};
    const std::vector<const methods::Method*> inexpensive_methods{&inexpensive_method};
    const auto inexpensive_result =
        methods::find_applicable_methods({.molecules = prepared_collection,
                                          .methods = inexpensive_methods,
                                          .parameter_sets = {},
                                          .resource_policy = {.cutoff_atom_threshold = 1}});
    REQUIRE(inexpensive_result.applicable.size() == 1);
    CHECK(assessment_for(inexpensive_result.applicable[0], calculation::ExecutionMode::full)
              .availability == methods::ExecutionAvailability::available);

    const chargefw::test::StubMethod expensive_method{
        "resource-test",
        {.resources = {.time = methods::ComplexityTerm::atoms_cubed,
                       .memory = methods::ComplexityTerm::atoms_squared}}};
    const std::vector<const methods::Method*> expensive_methods{&expensive_method};
    const auto below_threshold_result =
        methods::find_applicable_methods({.molecules = prepared_collection,
                                          .methods = expensive_methods,
                                          .parameter_sets = {},
                                          .resource_policy = {.cutoff_atom_threshold = 3}});
    REQUIRE(below_threshold_result.applicable.size() == 1);
    const auto& below_threshold_full =
        assessment_for(below_threshold_result.applicable[0], calculation::ExecutionMode::full);
    CHECK(below_threshold_full.availability == methods::ExecutionAvailability::available);
    CHECK(below_threshold_full.issues.empty());

    const auto threshold_result =
        methods::find_applicable_methods({.molecules = prepared_collection,
                                          .methods = expensive_methods,
                                          .parameter_sets = {},
                                          .resource_policy = {.cutoff_atom_threshold = 2}});
    REQUIRE(threshold_result.applicable.size() == 1);
    const auto& threshold_full =
        assessment_for(threshold_result.applicable[0], calculation::ExecutionMode::full);
    CHECK(threshold_full.availability == methods::ExecutionAvailability::available_with_warning);
    REQUIRE(threshold_full.issues.size() == 1);
    CHECK(threshold_full.issues[0].kind ==
          methods::ExecutionIssueKind::resource_threshold_exceeded);
    CHECK(threshold_full.issues[0].molecule_index == std::optional<std::size_t>{0});
    CHECK(assessment_for(threshold_result.applicable[0], calculation::ExecutionMode::cutoff)
              .availability == methods::ExecutionAvailability::unsupported);
    CHECK(assessment_for(threshold_result.applicable[0], calculation::ExecutionMode::cover)
              .availability == methods::ExecutionAvailability::unsupported);

    const std::vector<const methods::Method*> topology_methods{mgc};
    const auto topology_collection = chargefw::core::MoleculeCollection{
        std::vector{chargefw::test::make_water_graph()}, "topology-test-collection"};
    const features::PreparedMoleculeCollection topology_prepared{topology_collection};
    const auto topology_result =
        methods::find_applicable_methods({.molecules = topology_prepared,
                                          .methods = topology_methods,
                                          .parameter_sets = {},
                                          .resource_policy = {.cutoff_atom_threshold = 2}});
    REQUIRE(topology_result.applicable.size() == 1);
    CHECK(assessment_for(topology_result.applicable[0], calculation::ExecutionMode::full)
              .availability == methods::ExecutionAvailability::available_with_warning);
    CHECK(assessment_for(topology_result.applicable[0], calculation::ExecutionMode::cutoff)
              .availability == methods::ExecutionAvailability::unsupported);
    CHECK(assessment_for(topology_result.applicable[0], calculation::ExecutionMode::cover)
              .availability == methods::ExecutionAvailability::unsupported);

    const auto collection_threshold_result =
        methods::find_applicable_methods({.molecules = prepared_collection,
                                          .methods = expensive_methods,
                                          .parameter_sets = {},
                                          .resource_policy = {.cutoff_atom_threshold = 1}});
    REQUIRE(collection_threshold_result.applicable.size() == 1);
    const auto& collection_threshold_full =
        assessment_for(collection_threshold_result.applicable[0], calculation::ExecutionMode::full);
    CHECK(collection_threshold_full.issues.size() == prepared_collection.size());

    const auto unlimited_result = methods::find_applicable_methods(
        {.molecules = prepared_collection,
         .methods = expensive_methods,
         .parameter_sets = {},
         .resource_policy = {.cutoff_atom_threshold = std::nullopt,
                             .cover_atom_threshold = std::nullopt}});
    REQUIRE(unlimited_result.applicable.size() == 1);
    const auto& unlimited_full =
        assessment_for(unlimited_result.applicable[0], calculation::ExecutionMode::full);
    CHECK(unlimited_full.availability == methods::ExecutionAvailability::available);
    CHECK(unlimited_full.issues.empty());
}

TEST_CASE("reduced execution applicability requires geometry and declared support",
          "[methods][method-applicability]") {
    const auto collection = make_collection();
    const features::PreparedMoleculeCollection prepared_collection{collection};

    const chargefw::test::StubMethod topology_reduced_method{
        "resource-test", {.resources = {.supports_cutoff = true, .supports_cover = true}}};
    const std::vector<const methods::Method*> topology_reduced_methods{&topology_reduced_method};
    const auto topology_reduced_result =
        methods::find_applicable_methods({.molecules = prepared_collection,
                                          .methods = topology_reduced_methods,
                                          .parameter_sets = {}});
    REQUIRE(topology_reduced_result.applicable.size() == 1);
    CHECK(assessment_for(topology_reduced_result.applicable[0], calculation::ExecutionMode::cutoff)
              .availability == methods::ExecutionAvailability::unsupported);
    CHECK(assessment_for(topology_reduced_result.applicable[0], calculation::ExecutionMode::cover)
              .availability == methods::ExecutionAvailability::unsupported);

    const chargefw::test::StubMethod spatial_reduced_method{
        "resource-test",
        {.coordinates = true,
         .resources = {.supports_cutoff = true,
                       .supports_cover = true,
                       .reduced_charge_policy =
                           methods::ReducedChargePolicy::uniform_target_global}}};
    const std::vector<const methods::Method*> spatial_reduced_methods{&spatial_reduced_method};
    const core::MoleculeCollection spatial_collection{std::vector{chargefw::test::make_water()}};
    const features::PreparedMoleculeCollection spatial_prepared_collection{spatial_collection};
    const auto spatial_reduced_result =
        methods::find_applicable_methods({.molecules = spatial_prepared_collection,
                                          .methods = spatial_reduced_methods,
                                          .parameter_sets = {}});
    REQUIRE(spatial_reduced_result.applicable.size() == 1);
    CHECK(assessment_for(spatial_reduced_result.applicable[0], calculation::ExecutionMode::cutoff)
              .availability == methods::ExecutionAvailability::available);
    CHECK(assessment_for(spatial_reduced_result.applicable[0], calculation::ExecutionMode::cover)
              .availability == methods::ExecutionAvailability::available);
}
