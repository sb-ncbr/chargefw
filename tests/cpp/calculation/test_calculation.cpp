#include "support/test_calculation.h"
#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <algorithm>
#include <chargefw/calculation/calculation.h>
#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/core/position.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method.h>
#include <chargefw/methods/method_metadata.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/methods/method_registry.h>
#include <chargefw/methods/method_requirements.h>
#include <chargefw/parameters/models/common_parameters.h>
#include <chargefw/parameters/models/parameter_set.h>
#include <chargefw/parameters/models/parameter_set_metadata.h>
#include <cmath>
#include <limits>
#include <optional>
#include <snitch/snitch.hpp>
#include <span>
#include <stdexcept>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace calculation = chargefw::calculation;
namespace charges = chargefw::charges;
namespace core = chargefw::core;
namespace features = chargefw::features;
namespace methods = chargefw::methods;

using chargefw::test::calculate_application;

namespace {

auto make_parameter_set(std::string id, std::string method_id, const std::uint16_t priority)
    -> chargefw::parameters::ParameterSet {
    return chargefw::parameters::ParameterSet{
        chargefw::parameters::ParameterSetMetadata{.id = std::move(id),
                                                   .method_id = std::move(method_id),
                                                   .name = "Test parameters",
                                                   .priority = priority},
        {},
        chargefw::parameters::AtomParameters{
            {{.key = chargefw::test::atom_key(
                  1, chargefw::parameters::AtomParameterClassificationKind::PLAIN, "*"),
              .parameters = {{.name = "value", .value = 1.0}}},
             {.key = chargefw::test::atom_key(
                  8, chargefw::parameters::AtomParameterClassificationKind::PLAIN, "*"),
              .parameters = {{.name = "value", .value = 1.0}}}}}};
}

auto make_singular_eem_parameter_set() -> chargefw::parameters::ParameterSet {
    return chargefw::parameters::ParameterSet{
        chargefw::parameters::ParameterSetMetadata{
            .id = "singular-eem", .method_id = "eem", .name = "Singular EEM parameters"},
        chargefw::parameters::CommonParameters{{{.name = "kappa", .value = 1.0}}},
        chargefw::parameters::AtomParameters{
            {{.key = chargefw::test::atom_key(
                  1, chargefw::parameters::AtomParameterClassificationKind::PLAIN, "*"),
              .parameters = {{.name = "A", .value = 1.0}, {.name = "B", .value = 1.0}}},
             {.key = chargefw::test::atom_key(
                  8, chargefw::parameters::AtomParameterClassificationKind::PLAIN, "*"),
              .parameters = {{.name = "A", .value = 2.0}, {.name = "B", .value = 1.0}}}}}};
}

auto make_singular_eem_molecule() -> core::Molecule {
    return core::Molecule{
        std::vector{core::Atom{1}, core::Atom{8}},
        {},
        {core::Conformer{{core::Position{0.0, 0.0, 0.0}, core::Position{1.0, 0.0, 0.0}}}},
        "singular-eem"};
}

auto make_duplicate_parameter_request(const std::string_view first_method_id,
                                      const std::string_view second_method_id,
                                      const bool explicit_selection)
    -> calculation::AssessmentRequest {
    auto request = calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .parameter_sets = {make_parameter_set("duplicate", std::string{first_method_id}, 0),
                           make_parameter_set("duplicate", std::string{second_method_id}, 0)}};
    if (explicit_selection) {
        request.method_id = "qeq";
        request.parameter_set_id = "duplicate";
    }
    return request;
}

auto make_double_bonded_carbons() -> core::Molecule {
    return core::Molecule{std::vector{core::Atom{6}, core::Atom{6}},
                          std::vector{core::Bond{0, 1, core::BondOrder::DOUBLE}},
                          {},
                          "double-bonded-carbons"};
}

auto make_permissive_peoe_parameter_set() -> chargefw::parameters::ParameterSet {
    return chargefw::parameters::ParameterSet{
        chargefw::parameters::ParameterSetMetadata{.id = "permissive-peoe-parameters",
                                                   .method_id = "peoe",
                                                   .name = "Permissive PEOE parameters"},
        chargefw::parameters::CommonParameters{{{.name = "dampH", .value = 1.0}}},
        chargefw::parameters::AtomParameters{
            {{.key = chargefw::test::atom_key(
                  6, chargefw::parameters::AtomParameterClassificationKind::HIGHEST_BOND_ORDER,
                  "1"),
              .parameters = {{.name = "A", .value = 1.0},
                             {.name = "B", .value = 1.0},
                             {.name = "C", .value = 1.0}}}}}};
}

} // namespace

static_assert(std::is_move_constructible_v<calculation::AssessmentResult>);
static_assert(!std::is_move_assignable_v<calculation::AssessmentResult>);
static_assert(
    std::is_same_v<decltype(std::declval<const calculation::AssessmentResult&>().molecules()),
                   const core::MoleculeCollection&>);

TEST_CASE("assessment owns its inputs across lvalue, rvalue, and relocation",
          "[calculation][calculation]") {
    const auto make_request = [] {
        return calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::vector{make_double_bonded_carbons()}},
            .parameter_sets = {make_permissive_peoe_parameter_set()},
            .method_id = "peoe",
            .parameter_set_id = "permissive-peoe-parameters",
            .classification_options = {.permissive_types = true}};
    };
    // The lvalue overload copies its inputs, so the assessment outlives the request.
    const auto copied = [&] {
        const auto request = make_request();
        return calculation::assess(request);
    }();
    auto consumed = calculation::assess(make_request());
    const auto& owned_molecules = consumed.molecules();
    const auto relocated = std::move(consumed);
    CHECK(&relocated.molecules() == &owned_molecules);

    const auto copied_result = calculation::calculate(copied);
    const auto relocated_result = calculation::calculate(relocated);
    for (const auto* result : {&copied_result, &relocated_result}) {
        REQUIRE(result->calculated());
        REQUIRE(result->effective.has_value());
        CHECK(result->effective->method_id == "peoe");
        CHECK(result->effective->parameter_set_id ==
              std::optional<std::string>{"permissive-peoe-parameters"});
        CHECK(result->charges->assignment(0).charges.size() == 2);
    }
    CHECK(std::ranges::equal(copied_result.charges->assignment(0).charges.values(),
                             relocated_result.charges->assignment(0).charges.values()));
}

TEST_CASE("calculation facade reports singular solver failures with target context",
          "[calculation][calculation]") {
    const auto result = calculate_application(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{make_singular_eem_molecule()}},
        .parameter_sets = {make_singular_eem_parameter_set()},
        .method_id = "eem",
        .parameter_set_id = "singular-eem"});

    CHECK(result.status == calculation::ExecutionStatus::numerical_failure);
    CHECK(!result.charges.has_value());
    REQUIRE(result.effective.has_value());
    CHECK(result.effective->method_id == "eem");
    CHECK(result.effective->parameter_set_id == std::optional<std::string>{"singular-eem"});
    REQUIRE(result.failure_message.has_value());
    CHECK(
        result.failure_message->contains("method 'eem', molecule 1 ('singular-eem'), conformer 1"));
}

TEST_CASE("fixed-charge EEM retains provenance on numerical failure",
          "[calculation][calculation]") {
    const auto molecule =
        core::Molecule{std::vector{core::Atom{1, 0}, core::Atom{12, 2}, core::Atom{8, 0}},
                       {},
                       {core::Conformer{{{0.0, 0.0, 0.0}, {0.0, 3.0, 0.0}, {1.0, 0.0, 0.0}}}},
                       "singular-fixed-charge-eem"};
    const auto result = calculation::calculate(calculation::assess(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{molecule}},
        .parameter_sets = {make_singular_eem_parameter_set()},
        .method_id = "eem",
        .parameter_set_id = "singular-eem",
        .fixed_ions = calculation::FixedIons{
            .sources = {{.molecule_index = 0, .atom_index = 1, .charge = 0.25}}}}));

    CHECK(result.status == calculation::ExecutionStatus::numerical_failure);
    CHECK_FALSE(result.charges.has_value());
    REQUIRE(result.effective.has_value());
    REQUIRE(result.effective->fixed_ions.has_value());
    const auto& provenance = *result.effective->fixed_ions;
    REQUIRE(provenance.sources.size() == 1);
    CHECK(provenance.sources[0].charge == 0.25);
}

TEST_CASE("method options are validated and recorded in effective provenance",
          "[calculation][calculation]") {
    auto peoe_options = methods::MethodOptions{};
    peoe_options.set("iters", 1);
    const auto configured_peoe_result = calculate_application(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{make_double_bonded_carbons()}},
        .parameter_sets = {make_permissive_peoe_parameter_set()},
        .method_id = "peoe",
        .parameter_set_id = "permissive-peoe-parameters",
        .method_options = {{"peoe", peoe_options}},
        .classification_options = {.permissive_types = true}});
    REQUIRE(configured_peoe_result.calculated());
    REQUIRE(configured_peoe_result.effective.has_value());
    CHECK(configured_peoe_result.effective->method_options.get<int>("iters") == 1);

    auto invalid_peoe_options = methods::MethodOptions{};
    invalid_peoe_options.set("iters", std::string{"one"});
    const auto calculate_with_invalid_peoe_options = [&] -> void {
        static_cast<void>(calculate_application(calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::vector{make_double_bonded_carbons()}},
            .parameter_sets = {make_permissive_peoe_parameter_set()},
            .method_id = "peoe",
            .parameter_set_id = "permissive-peoe-parameters",
            .method_options = {{"peoe", invalid_peoe_options}},
            .classification_options = {.permissive_types = true}}));
    };
    CHECK_THROWS_AS(calculate_with_invalid_peoe_options(), std::invalid_argument);
}

TEST_CASE("assessment reports parameter and automatic method rejections",
          "[calculation][calculation]") {
    const auto rejected_parameter_assessment = calculation::assess(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{make_double_bonded_carbons()}},
        .parameter_sets = {make_permissive_peoe_parameter_set()},
        .method_id = "peoe",
        .parameter_set_id = "permissive-peoe-parameters"});
    CHECK(rejected_parameter_assessment.plans().empty());
    REQUIRE(rejected_parameter_assessment.rejections().size() == 1);
    CHECK(rejected_parameter_assessment.rejections()[0].method_id == "peoe");
    CHECK(rejected_parameter_assessment.rejections()[0].parameter_set_id ==
          std::optional<std::string>{"permissive-peoe-parameters"});

    auto automatic_rejected_assessment = calculation::assess(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .method_scope = calculation::MethodScope::all});
    const auto automatic_rejected_result = calculation::calculate(automatic_rejected_assessment);
    const auto rejected_smpqeq = std::ranges::find_if(
        automatic_rejected_result.rejections,
        [](const calculation::Rejection& rejection) { return rejection.method_id == "smpqeq"; });
    REQUIRE(rejected_smpqeq != automatic_rejected_result.rejections.end());
    CHECK(!rejected_smpqeq->issues.empty());
}

TEST_CASE("automatic selection considers methods that are not recommended only when asked for all",
          "[calculation][calculation]") {
    const auto plan_method_ids = [](const calculation::AssessmentRequest& request) {
        const auto assessment = calculation::assess(request);
        auto ids = std::vector<std::string_view>{};
        for (const auto& plan : assessment.plans()) {
            ids.push_back(plan.method().id());
        }
        return ids;
    };
    const auto water = core::MoleculeCollection{std::vector{chargefw::test::make_water()}};

    // VEEM needs no parameters and applies to water, but it is not recommended.
    const auto recommended = plan_method_ids({.molecules = water});
    REQUIRE_FALSE(recommended.empty());
    CHECK(std::ranges::find(recommended, "veem") == recommended.end());

    const auto all =
        plan_method_ids({.molecules = water, .method_scope = calculation::MethodScope::all});
    CHECK(std::ranges::find(all, "veem") != all.end());

    // A named method is assessed whether or not it is recommended.
    const auto named = plan_method_ids({.molecules = water, .method_id = "veem"});
    REQUIRE(named.size() == 1);
    CHECK(named[0] == "veem");
}

TEST_CASE("explicit no-plan assessment reports rejected scientific prerequisites",
          "[calculation][calculation]") {
    auto assessment = calculation::assess(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .method_id = "smpqeq"});

    CHECK(assessment.plans().empty());
    REQUIRE(assessment.rejections().size() == 1);
    CHECK_FALSE(assessment.rejections()[0].policy.has_value());
    CHECK(assessment.rejections()[0].method_id == "smpqeq");
    CHECK_FALSE(assessment.rejections()[0].issues.empty());

    const auto result = calculation::calculate(assessment);
    CHECK(result.status == calculation::ExecutionStatus::no_executable_plan);
    CHECK_FALSE(result.calculated());
}

TEST_CASE("explicit full execution reports its plan and effective provenance",
          "[calculation][calculation]") {
    const auto assessment = calculation::assess(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .method_id = "formal",
        .execution_selection =
            calculation::ExecutionSelection{calculation::ExecutionSelectionKind::full}});
    REQUIRE(assessment.plans().size() == 1);
    REQUIRE(assessment.default_plan() != nullptr);
    CHECK(assessment.default_plan()->method().id() == "formal");
    CHECK(assessment.default_plan()->policy().mode() == calculation::ExecutionMode::full);

    const auto result = calculation::calculate(assessment, 1);
    REQUIRE(result.calculated());
    CHECK(result.charges->method_id() == std::string_view{"formal"});
    CHECK(result.charges->size() == 1);
    REQUIRE(result.effective.has_value());
    CHECK(result.effective->method_id == "formal");
    CHECK_FALSE(result.effective->parameter_set_id.has_value());
    CHECK(result.effective->execution_policy.mode() == calculation::ExecutionMode::full);
    CHECK(result.effective->execution_issues.empty());
    CHECK(result.metrics.applicability_seconds >= 0.0);
}

TEST_CASE("resource thresholds guide automatic execution and warn on explicit full execution",
          "[calculation][calculation]") {
    const auto automatic_fallback_result = calculate_application(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .parameter_sets = {},
        .resource_policy = {.cutoff_atom_threshold = 2}});
    REQUIRE(automatic_fallback_result.calculated());
    CHECK(automatic_fallback_result.charges->method_id() == std::string_view{"eqeq"});
    REQUIRE(automatic_fallback_result.effective.has_value());
    CHECK(automatic_fallback_result.effective->method_id == "eqeq");
    CHECK(automatic_fallback_result.effective->execution_policy.mode() ==
          calculation::ExecutionMode::cutoff);
    CHECK(automatic_fallback_result.effective->execution_issues.empty());

    const auto automatic_mgc_result = calculate_application(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .parameter_sets = {},
        .method_id = "mgc",
        .resource_policy = {.cutoff_atom_threshold = 2}});
    CHECK(automatic_mgc_result.status == calculation::ExecutionStatus::no_executable_plan);
    CHECK_FALSE(automatic_mgc_result.calculated());
    CHECK_FALSE(automatic_mgc_result.effective.has_value());

    const auto explicit_full_result = calculate_application(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .parameter_sets = {},
        .method_id = "mgc",
        .execution_selection =
            calculation::ExecutionSelection{calculation::ExecutionSelectionKind::full},
        .resource_policy = {.cutoff_atom_threshold = 2}});
    REQUIRE(explicit_full_result.calculated());
    CHECK(explicit_full_result.charges->method_id() == std::string_view{"mgc"});
    REQUIRE(explicit_full_result.effective.has_value());
    CHECK(explicit_full_result.effective->execution_policy.mode() ==
          calculation::ExecutionMode::full);
    REQUIRE(explicit_full_result.effective->execution_issues.size() == 1);
    CHECK(explicit_full_result.effective->execution_issues[0].kind ==
          methods::ExecutionIssueKind::resource_threshold_exceeded);

    const auto unlimited_result = calculate_application(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .parameter_sets = {},
        .method_id = "mgc",
        .resource_policy = {.cutoff_atom_threshold = std::nullopt,
                            .cover_atom_threshold = std::nullopt}});
    REQUIRE(unlimited_result.calculated());
    CHECK(unlimited_result.charges->method_id() == std::string_view{"mgc"});
    REQUIRE(unlimited_result.effective.has_value());
    CHECK(unlimited_result.effective->execution_issues.empty());

    const auto unsupported_cutoff_result = calculate_application(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .parameter_sets = {},
        .method_id = "mgc",
        .execution_selection =
            calculation::ExecutionSelection{calculation::ExecutionSelectionKind::cutoff, 8.0}});
    CHECK(unsupported_cutoff_result.status == calculation::ExecutionStatus::no_executable_plan);
    CHECK_FALSE(unsupported_cutoff_result.calculated());
}

TEST_CASE("plan estimates compose unit costs by execution mode", "[calculation][calculation]") {
    // Four waters 20 A apart: every cutoff fragment and cover pivot holds exactly one water.
    auto atoms = std::vector<core::Atom>{};
    auto bonds = std::vector<core::Bond>{};
    auto positions = std::vector<core::Position>{};
    const auto water = chargefw::test::make_water();
    for (std::size_t copy = 0; copy < 4; ++copy) {
        const auto offset = atoms.size();
        for (std::size_t atom_index = 0; atom_index < water.atom_count(); ++atom_index) {
            atoms.push_back(water.atom(atom_index));
            auto position = water.conformer(0).positions()[atom_index];
            position.x += 20.0 * static_cast<double>(copy);
            positions.push_back(position);
        }
        for (const auto& bond : water.bonds()) {
            bonds.emplace_back(offset + bond.first_atom_index(), offset + bond.second_atom_index());
        }
    }
    const auto cluster = core::Molecule{
        std::move(atoms), std::move(bonds), {core::Conformer{std::move(positions)}}, "waters"};
    const auto molecules = core::MoleculeCollection{std::vector{water, cluster}};
    const auto* eqeq = methods::method_registry().find("eqeq");
    REQUIRE(eqeq != nullptr);
    const auto options = methods::make_default_options(eqeq->option_schema());
    const auto water_cost = eqeq->estimate({.atoms = 3, .bonds = 2}, options);
    const auto cluster_cost = eqeq->estimate({.atoms = 12, .bonds = 8}, options);

    const auto full = calculation::assess(
        calculation::AssessmentRequest{.molecules = molecules, .method_id = "eqeq"});
    REQUIRE(full.default_plan() != nullptr);
    REQUIRE(full.default_plan()->policy().mode() == calculation::ExecutionMode::full);
    const auto serial = calculation::estimate(full, *full.default_plan(), 1);
    const auto parallel = calculation::estimate(full, *full.default_plan(), 2);
    const auto oversubscribed = calculation::estimate(full, *full.default_plan(), 8);
    CHECK(std::abs(serial.operations - (water_cost.operations + cluster_cost.operations)) <
          1e-9 * serial.operations);
    // A second thread adds the smaller target's workspace; further threads have no targets to run.
    CHECK(std::abs(parallel.memory_bytes - serial.memory_bytes - water_cost.memory_bytes) <
          1e-9 * parallel.memory_bytes);
    CHECK(oversubscribed.memory_bytes == parallel.memory_bytes);
    CHECK(calculation::estimate(full, *full.default_plan()).memory_bytes >= serial.memory_bytes);
    // A limit above oneTBB's allowed parallelism runs and estimates at that parallelism.
    const auto unlimited = std::numeric_limits<std::size_t>::max();
    CHECK(calculation::estimate(full, *full.default_plan(), unlimited).memory_bytes ==
          calculation::estimate(full, *full.default_plan()).memory_bytes);
    CHECK(calculation::calculate(full, *full.default_plan(), unlimited).calculated());

    const auto reduced = [&](const calculation::ExecutionSelectionKind kind) {
        return calculation::assess(calculation::AssessmentRequest{
            .molecules = molecules,
            .method_id = "eqeq",
            .execution_selection =
                calculation::ExecutionSelection{kind, calculation::minimum_reduced_radius}});
    };
    const auto cutoff = reduced(calculation::ExecutionSelectionKind::cutoff);
    const auto cover = reduced(calculation::ExecutionSelectionKind::cover);
    REQUIRE(cutoff.default_plan() != nullptr);
    REQUIRE(cover.default_plan() != nullptr);
    const auto cutoff_serial = calculation::estimate(cutoff, *cutoff.default_plan(), 1);
    const auto cutoff_parallel = calculation::estimate(cutoff, *cutoff.default_plan(), 4);
    const auto cover_serial = calculation::estimate(cover, *cover.default_plan(), 1);
    // Cutoff runs one fragment per atom in parallel; cover runs fewer fragments.
    CHECK(cutoff_parallel.memory_bytes > cutoff_serial.memory_bytes);
    CHECK(cover_serial.operations < cutoff_serial.operations);
}

TEST_CASE("assessment rejects invalid requests", "[calculation][calculation]") {
    const auto water = core::MoleculeCollection{std::vector{chargefw::test::make_water()}};
    const auto invalid_requests = std::vector<calculation::AssessmentRequest>{
        {.molecules = water, .method_id = "missing"},
        {.molecules = water, .method_id = "formal", .parameter_set_id = "missing"},
        {.molecules = water,
         .parameter_sets = {make_singular_eem_parameter_set()},
         .parameter_set_id = "singular-eem"},
        {.molecules = water,
         .parameter_sets = {make_singular_eem_parameter_set()},
         .method_id = "formal",
         .parameter_set_id = "singular-eem"},
        {.molecules = water,
         .parameter_sets = {make_singular_eem_parameter_set()},
         .method_id = "qeq",
         .parameter_set_id = "singular-eem"},
        {.molecules = water,
         .resource_policy = {.cutoff_atom_threshold = std::nullopt, .cover_atom_threshold = 10}},
        {.molecules = water,
         .resource_policy = {.cutoff_atom_threshold = 20, .cover_atom_threshold = 10}}};
    for (const auto& request : invalid_requests) {
        CAPTURE(request.method_id.value_or("automatic"), request.parameter_set_id.value_or(""));
        CHECK_THROWS_AS(static_cast<void>(calculation::assess(request)), std::invalid_argument);
    }
}

TEST_CASE("assessments expose reusable target-bound execution plans",
          "[calculation][calculation]") {
    auto assessment = calculation::assess(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .method_id = "eqeq"});

    REQUIRE_FALSE(assessment.plans().empty());
    REQUIRE(assessment.default_plan() == &assessment.plans().front());
    CHECK(assessment.plans().front().method().id() == "eqeq");
    CHECK(assessment.plans().front().policy().mode() == calculation::ExecutionMode::full);
    CHECK(assessment.plans().size() >= 2);

    const auto first = calculation::calculate(assessment, assessment.plans().front(), 1);
    const auto repeated = calculation::calculate(assessment, assessment.plans().front(), 1);
    REQUIRE(first.calculated());
    REQUIRE(repeated.calculated());
    CHECK(std::ranges::equal(first.charges->assignment(0).charges.values(),
                             repeated.charges->assignment(0).charges.values()));

    const auto cutoff_only = calculation::assess(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .method_id = "eqeq",
        .execution_selection = calculation::ExecutionSelection{
            calculation::ExecutionSelectionKind::cutoff, calculation::minimum_reduced_radius}});
    REQUIRE_FALSE(cutoff_only.plans().empty());
    CHECK(std::ranges::all_of(cutoff_only.plans(), [](const calculation::ExecutionPlan& plan) {
        return plan.policy().mode() == calculation::ExecutionMode::cutoff;
    }));

    const auto resource_limited = calculation::assess(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .method_id = "eqeq",
        .resource_policy = {.cutoff_atom_threshold = 2}});
    REQUIRE_FALSE(resource_limited.plans().empty());
    CHECK(resource_limited.default_plan()->policy().mode() == calculation::ExecutionMode::cutoff);
    CHECK(
        std::ranges::none_of(resource_limited.plans(), [](const calculation::ExecutionPlan& plan) {
            return plan.policy().mode() == calculation::ExecutionMode::full;
        }));
    CHECK(std::ranges::any_of(
        resource_limited.rejections(), [](const calculation::Rejection& rejection) {
            return rejection.policy.has_value() &&
                   rejection.policy->mode() == calculation::ExecutionMode::full;
        }));

    auto other_assessment = calculation::assess(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .method_id = "eqeq"});
    const auto calculate_mismatched_plan = [&] {
        static_cast<void>(calculation::calculate(other_assessment, assessment.plans().front()));
    };
    CHECK_THROWS_AS(calculate_mismatched_plan(), std::invalid_argument);
}

TEST_CASE("calculation preserves empty-input cardinality", "[calculation][calculation]") {

    // Empty collections and empty targets retain their distinct cardinalities: no collection
    // entries produce no assignments, while an empty molecule still produces one source assignment.
    const auto empty_collection_result = calculate_application(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector<core::Molecule>{}},
        .method_id = "formal"});
    REQUIRE(empty_collection_result.calculated());
    CHECK(empty_collection_result.charges->empty());

    const auto empty_target_result = calculate_application(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{core::Molecule{{}, {}, {}, "empty"}}},
        .method_id = "formal"});
    REQUIRE(empty_target_result.calculated());
    REQUIRE(empty_target_result.charges->size() == 1);
    CHECK(empty_target_result.charges->assignment(0).charges.empty());
}

TEST_CASE("both assessment overloads reject duplicate parameter-set IDs",
          "[calculation][calculation]") {
    // AssessmentRequest rejects duplicate parameter-set IDs before filtering or applicability. This
    // remains true when the duplicates target the same or different methods and for both overloads.
    for (const auto explicit_selection : {false, true}) {
        for (const auto& [first_method_id, second_method_id] :
             {std::pair{"qeq", "qeq"}, std::pair{"qeq", "eem"}}) {
            const auto lvalue_request = make_duplicate_parameter_request(
                first_method_id, second_method_id, explicit_selection);
            const auto assess_duplicate_lvalue_request = [&lvalue_request] -> void {
                static_cast<void>(calculation::assess(lvalue_request));
            };
            CHECK_THROWS_AS(assess_duplicate_lvalue_request(), std::invalid_argument);

            auto rvalue_request = make_duplicate_parameter_request(
                first_method_id, second_method_id, explicit_selection);
            const auto assess_duplicate_rvalue_request = [&rvalue_request] -> void {
                static_cast<void>(calculation::assess(std::move(rvalue_request)));
            };
            CHECK_THROWS_AS(assess_duplicate_rvalue_request(), std::invalid_argument);
        }
    }
}
