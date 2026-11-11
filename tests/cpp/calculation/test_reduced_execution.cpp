#include "calculation/calculation_request.h"
#include "calculation/cover_execution.h"
#include "calculation/cutoff_execution.h"
#include "calculation/fixed_charge_partition.h"
#include "calculation/reduced_execution.h"
#include "features/topology_helpers.h"
#include "support/test_calculation.h"
#include "support/test_methods.h"
#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <chargefw/calculation/calculation.h>
#include <chargefw/calculation/observer.h>
#include <chargefw/charges/atomic_charges.h>
#include <chargefw/core/atom.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/core/position.h>
#include <chargefw/features/conformer_features.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/features/spatial_fragment.h>
#include <chargefw/methods/calculation_input.h>
#include <chargefw/methods/method.h>
#include <chargefw/methods/method_metadata.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/methods/method_registry.h>
#include <chargefw/methods/method_requirements.h>
#include <chargefw/parameters/classification/parameter_classifier.h>
#include <chargefw/parameters/models/atom_parameters.h>
#include <chargefw/parameters/models/bond_parameters.h>
#include <chargefw/parameters/models/common_parameters.h>
#include <chargefw/parameters/models/parameter_key.h>
#include <chargefw/parameters/models/parameter_set.h>
#include <chargefw/parameters/models/parameter_set_metadata.h>
#include <chargefw/parameters/models/parameter_view.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <optional>
#include <snitch/snitch.hpp>
#include <span>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

namespace calculation = chargefw::calculation;
namespace charges = chargefw::charges;
namespace core = chargefw::core;
namespace features = chargefw::features;
namespace methods = chargefw::methods;
namespace parameters = chargefw::parameters;

using chargefw::test::calculate_application;

namespace {

class FailingFixedIonsFragmentMethod final : public chargefw::test::StubMethod {
  public:
    FailingFixedIonsFragmentMethod()
        : StubMethod{"failing-fixed-ions-fragment",
                     {.coordinates = true,
                      .resources = {.supports_cutoff = true,
                                    .supports_cover = true,
                                    .reduced_charge_policy =
                                        methods::ReducedChargePolicy::uniform_target_global},
                      .supports_fixed_point_sources = true}} {}

    [[nodiscard]] auto calculate(const methods::CalculationInput&) const
        -> charges::AtomicCharges override {
        throw std::runtime_error{"deliberate fragment failure"};
    }
};

auto make_eqeqc_parameters() -> parameters::ParameterSet {
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{
            .id = "test-eqeqc", .method_id = "eqeqc", .name = "Test EQeq+C"},
        parameters::CommonParameters{{{.name = "alpha", .value = 1.0}}},
        parameters::AtomParameters{{{.key = chargefw::test::plain_atom_key(1),
                                     .parameters = {{.name = "Dz", .value = 0.1}}},
                                    {.key = chargefw::test::plain_atom_key(8),
                                     .parameters = {{.name = "Dz", .value = 0.2}}}}}};
}

auto make_charged_water() -> core::Molecule {
    return core::Molecule{
        std::vector{core::Atom{8, 1, "O"}, core::Atom{1, 0, "H1"}, core::Atom{1, 0, "H2"}},
        std::vector{core::Bond{0, 1, core::BondOrder::SINGLE},
                    core::Bond{0, 2, core::BondOrder::SINGLE}},
        std::vector{core::Conformer{{core::Position{.x = 0.0, .y = 0.0, .z = 0.0},
                                     core::Position{.x = 0.9572, .y = 0.0, .z = 0.0},
                                     core::Position{.x = -0.2390, .y = 0.9270, .z = 0.0}}}},
        "charged-water"};
}

auto make_extended_components(const int isolated_formal_charge) -> core::Molecule {
    constexpr auto chain_atom_count = std::size_t{12};
    auto atoms = std::vector<core::Atom>{};
    auto bonds = std::vector<core::Bond>{};
    auto positions = std::vector<core::Position>{};
    atoms.reserve(chain_atom_count + 1);
    bonds.reserve(chain_atom_count - 1);
    positions.reserve(chain_atom_count + 1);

    for (std::size_t atom_index = 0; atom_index < chain_atom_count; ++atom_index) {
        const auto atomic_number = atom_index % 2 == 0 ? 8 : 1;
        atoms.emplace_back(atomic_number, 0, "chain");
        positions.push_back(core::Position{.x = 1.4 * static_cast<double>(atom_index)});
        if (atom_index != 0) {
            bonds.emplace_back(atom_index - 1, atom_index, core::BondOrder::SINGLE);
        }
    }
    atoms.emplace_back(1, isolated_formal_charge, "isolated");
    positions.push_back(core::Position{.x = 30.0});

    return core::Molecule{std::move(atoms), std::move(bonds),
                          std::vector{core::Conformer{std::move(positions)}},
                          "extended-components"};
}

auto make_two_diatomic_components(const double separation) -> core::Molecule {
    return core::Molecule{
        std::vector{core::Atom{8, 0, "O1"}, core::Atom{1, 0, "H1"}, core::Atom{8, 0, "O2"},
                    core::Atom{1, 0, "H2"}},
        std::vector{core::Bond{0, 1, core::BondOrder::SINGLE},
                    core::Bond{2, 3, core::BondOrder::SINGLE}},
        std::vector{core::Conformer{{core::Position{.x = 0.0}, core::Position{.x = 1.0},
                                     core::Position{.x = separation},
                                     core::Position{.x = separation + 1.2}}}},
        "two-diatomic-components"};
}

[[nodiscard]] auto calculate_reduced(const core::Molecule& molecule,
                                     const std::string_view method_id,
                                     std::vector<parameters::ParameterSet> parameter_sets,
                                     const calculation::ExecutionSelectionKind mode,
                                     const double radius = 8.0) -> charges::AtomicCharges {
    const auto result = calculate_application(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{molecule}},
        .parameter_sets = std::move(parameter_sets),
        .method_id = std::string{method_id},
        .execution_selection = calculation::ExecutionSelection{mode, radius}});
    REQUIRE(result.calculated());
    REQUIRE(result.charges.has_value());
    return result.charges->assignment(0).charges;
}

[[nodiscard]] auto calculate_full(const core::Molecule& molecule, const std::string_view method_id,
                                  std::vector<parameters::ParameterSet> parameter_sets)
    -> charges::AtomicCharges {
    const auto result = calculate_application(
        calculation::AssessmentRequest{.molecules = core::MoleculeCollection{std::vector{molecule}},
                                       .parameter_sets = std::move(parameter_sets),
                                       .method_id = std::string{method_id}});
    REQUIRE(result.calculated());
    REQUIRE(result.charges.has_value());
    return result.charges->assignment(0).charges;
}

struct ErrorMetrics {
    double rmsd;
    double mae;
    double maxabs;
};

[[nodiscard]] auto error_metrics(const charges::AtomicCharges& actual,
                                 const charges::AtomicCharges& reference) -> ErrorMetrics {
    REQUIRE(actual.size() == reference.size());
    auto squared_error = 0.0;
    auto absolute_error = 0.0;
    auto maxabs = 0.0;
    for (std::size_t atom_index = 0; atom_index < actual.size(); ++atom_index) {
        const auto difference = actual[atom_index] - reference[atom_index];
        const auto absolute_difference = std::abs(difference);
        squared_error += difference * difference;
        absolute_error += absolute_difference;
        maxabs = std::max(maxabs, absolute_difference);
    }
    const auto atom_count = static_cast<double>(actual.size());
    return {.rmsd = std::sqrt(squared_error / atom_count),
            .mae = absolute_error / atom_count,
            .maxabs = maxabs};
}

auto check_component_total(const charges::AtomicCharges& atomic_charges, const std::size_t begin,
                           const std::size_t end, const double expected) -> void {
    auto total = 0.0;
    for (auto atom_index = begin; atom_index < end; ++atom_index) {
        total += atomic_charges[atom_index];
    }
    CHECK(std::abs(total - expected) < 1.0e-10);
}

auto assert_reduced_matches_full(
    const std::string_view method_id,
    const std::vector<parameters::ParameterSet>& parameter_sets = {},
    core::Molecule molecule = chargefw::test::make_two_conformer_water()) -> void {
    const auto molecules = core::MoleculeCollection{std::vector{molecule, std::move(molecule)}};
    const auto full =
        calculate_application(calculation::AssessmentRequest{.molecules = molecules,
                                                             .parameter_sets = parameter_sets,
                                                             .method_id = std::string{method_id}});
    REQUIRE(full.calculated());

    for (const auto selection_kind : {calculation::ExecutionSelectionKind::cutoff,
                                      calculation::ExecutionSelectionKind::cover}) {
        const auto reduced = calculate_application(calculation::AssessmentRequest{
            .molecules = molecules,
            .parameter_sets = parameter_sets,
            .method_id = std::string{method_id},
            .execution_selection = calculation::ExecutionSelection{selection_kind, 8.0}});

        REQUIRE(reduced.calculated());
        REQUIRE(reduced.effective.has_value());
        CHECK(reduced.effective->execution_policy.mode() ==
              (selection_kind == calculation::ExecutionSelectionKind::cutoff
                   ? calculation::ExecutionMode::cutoff
                   : calculation::ExecutionMode::cover));
        REQUIRE(full.charges->size() == reduced.charges->size());

        for (std::size_t assignment_index = 0; assignment_index < full.charges->size();
             ++assignment_index) {
            const auto& full_charges = full.charges->assignment(assignment_index).charges;
            const auto& reduced_charges = reduced.charges->assignment(assignment_index).charges;
            chargefw::test::assert_same_charges(reduced_charges, full_charges, 1.0e-10);
        }
    }
}

} // namespace

TEST_CASE("fragment classification projects source entries to local indices",
          "[calculation][reduced]") {
    const auto molecule = make_two_diatomic_components(10.0);
    const features::PreparedMolecule prepared{molecule};
    const features::ConformerFeatures geometry{molecule, 0};
    const auto fragment = features::SpatialFragmentBuilder{prepared, geometry}.build(2, 1.2);
    REQUIRE(std::ranges::equal(fragment.local_to_source_atom_indices(),
                               std::vector<std::size_t>{2, 3}));

    const auto projected = calculation::detail::project_classification(
        parameters::ParameterClassification{
            parameters::AtomParameterClassification{{10, 11, 12, 13}},
            parameters::BondParameterClassification{{20, 21}}},
        fragment);
    CHECK(std::ranges::equal(projected.atom().parameter_entry_indices(),
                             std::vector<std::size_t>{12, 13}));
    CHECK(std::ranges::equal(projected.bond().parameter_entry_indices(),
                             std::vector<std::size_t>{21}));

    const auto atom_only = calculation::detail::project_classification(
        parameters::ParameterClassification{
            parameters::AtomParameterClassification{{10, 11, 12, 13}}},
        fragment);
    CHECK(std::ranges::equal(atom_only.atom().parameter_entry_indices(),
                             std::vector<std::size_t>{12, 13}));
    CHECK(atom_only.bond().empty());
}

TEST_CASE("direct reduced execution validates inputs and distributes a uniform target",
          "[calculation][reduced-execution]") {
    const chargefw::test::StubMethod zero_method{
        "zero-fragment",
        {.coordinates = true,
         .resources = {.supports_cutoff = true,
                       .supports_cover = true,
                       .reduced_charge_policy =
                           methods::ReducedChargePolicy::uniform_target_global}}};
    const auto charged_molecule = core::Molecule{
        std::vector{core::Atom{1, 1}, core::Atom{1, 0}},
        {},
        std::vector{core::Conformer{{core::Position{.x = 0.0}, core::Position{.x = 20.0}}}},
        "charged"};
    const auto collection = core::MoleculeCollection{std::vector{charged_molecule}};
    const features::PreparedMoleculeCollection prepared{collection};
    const methods::ApplicableMethod selected{.method = &zero_method, .parameter_set = nullptr};

    const methods::ApplicableMethod invalid_selected{
        .method = &zero_method, .parameter_set = nullptr, .classifications = {{}}};
    const auto no_conformer_collection = core::MoleculeCollection{std::vector{
        charged_molecule, core::Molecule{std::vector{core::Atom{1}}, {}, {}, "no-conformer"}}};
    const features::PreparedMoleculeCollection no_conformer_prepared{no_conformer_collection};
    for (const auto mode : {calculation::ExecutionMode::full, calculation::ExecutionMode::cutoff,
                            calculation::ExecutionMode::cover}) {
        CAPTURE(calculation::to_string(mode));
        const auto policy = mode == calculation::ExecutionMode::full
                                ? calculation::ExecutionPolicy{}
                                : calculation::ExecutionPolicy{mode, 8.0};
        CHECK_THROWS_AS(
            calculation::calculate(
                {.molecules = prepared, .selected = invalid_selected, .execution_policy = policy}),
            std::invalid_argument);
        CHECK_THROWS_MATCHES(calculation::calculate({.molecules = no_conformer_prepared,
                                                     .selected = selected,
                                                     .execution_policy = policy}),
                             std::invalid_argument,
                             snitch::matchers::with_what_contains{"molecule 2"});
    }

    const auto cutoff = calculation::calculate(
        {.molecules = prepared,
         .selected = selected,
         .execution_policy = calculation::ExecutionPolicy{calculation::ExecutionMode::cutoff, 8.0},
         .max_threads = 2});
    CHECK(cutoff.charges.assignment(0).charges[0] == 0.5);
    CHECK(cutoff.charges.assignment(0).charges[1] == 0.5);

    const auto cover = calculation::calculate(
        {.molecules = prepared,
         .selected = selected,
         .execution_policy = calculation::ExecutionPolicy{calculation::ExecutionMode::cover, 8.0}});
    CHECK(cover.charges.assignment(0).charges[0] == 0.5);
    CHECK(cover.charges.assignment(0).charges[1] == 0.5);
}

TEST_CASE("automatic execution selects reduced modes from resource thresholds",
          "[calculation][reduced-execution]") {
    const auto automatic_cutoff = calculate_application(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .parameter_sets = {chargefw::test::make_eem_parameters()},
        .method_id = "eem",
        .resource_policy = {.cutoff_atom_threshold = 2}});
    REQUIRE(automatic_cutoff.calculated());
    REQUIRE(automatic_cutoff.effective.has_value());
    CHECK(automatic_cutoff.effective->execution_policy.mode() ==
          calculation::ExecutionMode::cutoff);
    CHECK(automatic_cutoff.effective->execution_policy.radius() ==
          std::optional<double>{calculation::default_automatic_reduced_radius});

    const auto overridden_automatic_cutoff = calculate_application(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .parameter_sets = {chargefw::test::make_eem_parameters()},
        .method_id = "eem",
        .execution_selection =
            calculation::ExecutionSelection{calculation::ExecutionSelectionKind::automatic, 8.0},
        .resource_policy = {.cutoff_atom_threshold = 2}});
    REQUIRE(overridden_automatic_cutoff.calculated());
    REQUIRE(overridden_automatic_cutoff.effective.has_value());
    CHECK(overridden_automatic_cutoff.effective->execution_policy.mode() ==
          calculation::ExecutionMode::cutoff);
    CHECK(overridden_automatic_cutoff.effective->execution_policy.radius() ==
          std::optional<double>{8.0});

    const auto automatic_cover = calculate_application(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .parameter_sets = {chargefw::test::make_eem_parameters()},
        .method_id = "eem",
        .resource_policy = {.cutoff_atom_threshold = 2, .cover_atom_threshold = 2}});
    REQUIRE(automatic_cover.calculated());
    REQUIRE(automatic_cover.effective.has_value());
    CHECK(automatic_cover.effective->execution_policy.mode() == calculation::ExecutionMode::cover);

    const auto explicit_cutoff_above_cover_threshold =
        calculate_application(calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
            .parameter_sets = {chargefw::test::make_eem_parameters()},
            .method_id = "eem",
            .execution_selection =
                calculation::ExecutionSelection{calculation::ExecutionSelectionKind::cutoff,
                                                calculation::minimum_reduced_radius},
            .resource_policy = {.cutoff_atom_threshold = 2, .cover_atom_threshold = 2}});
    REQUIRE(explicit_cutoff_above_cover_threshold.calculated());
    REQUIRE(explicit_cutoff_above_cover_threshold.effective.has_value());
    CHECK(explicit_cutoff_above_cover_threshold.effective->execution_policy.mode() ==
          calculation::ExecutionMode::cutoff);
    CHECK(explicit_cutoff_above_cover_threshold.effective->execution_issues.size() == 1);
}

TEST_CASE("whole-molecule reduced execution matches full execution for each charge policy",
          "[calculation][reduced-execution]") {
    struct PolicyCase {
        methods::ReducedChargePolicy policy;
        std::string_view method_id;
        parameters::ParameterSet parameter_set;
        core::Molecule molecule;
    };
    using Policy = methods::ReducedChargePolicy;
    const auto water = chargefw::test::make_two_conformer_water();
    const auto point_charges =
        chargefw::test::SqeHoValues{.hydrogen_width = 0.0, .oxygen_width = 0.0};
    const auto reference_charges =
        chargefw::test::SqeHoValues{.hydrogen_q0 = 0.25, .oxygen_q0 = -0.5};
    const auto point_reference_charges = chargefw::test::SqeHoValues{
        .hydrogen_width = 0.0, .oxygen_width = 0.0, .hydrogen_q0 = 0.25, .oxygen_q0 = -0.5};
    const auto cases = std::vector<PolicyCase>{
        {Policy::uniform_target_global, "eem", chargefw::test::make_eem_parameters(), water},
        {Policy::zero_components, "sqe", chargefw::test::make_sqe_ho_parameters("sqe"), water},
        {Policy::zero_components, "sqe",
         chargefw::test::make_sqe_ho_parameters("sqe", point_charges), water},
        {Policy::formal_charge_components, "sqeq0", chargefw::test::make_sqe_ho_parameters("sqeq0"),
         water},
        {Policy::formal_charge_components, "sqeq0",
         chargefw::test::make_sqe_ho_parameters("sqeq0", point_charges), water},
        {Policy::formal_charge_components, "sqeq0", chargefw::test::make_sqe_ho_parameters("sqeq0"),
         make_charged_water()},
        {Policy::parameterized_charge_components, "sqeqp",
         chargefw::test::make_sqe_ho_parameters("sqeqp", reference_charges), water},
        {Policy::parameterized_charge_components, "sqeqp",
         chargefw::test::make_sqe_ho_parameters("sqeqp", point_reference_charges), water},
        {Policy::parameterized_charge_components, "sqeqp",
         chargefw::test::make_sqe_ho_parameters("sqeqp", reference_charges), make_charged_water()}};

    for (const auto& test_case : cases) {
        CAPTURE(test_case.method_id, test_case.molecule.name());
        const auto* method = methods::method_registry().find(test_case.method_id);
        REQUIRE(method != nullptr);
        CHECK(method->requirements().resources.reduced_charge_policy == test_case.policy);
        assert_reduced_matches_full(test_case.method_id, {test_case.parameter_set},
                                    test_case.molecule);
    }
}

TEST_CASE("SQE reduced execution preserves original component charge budgets",
          "[calculation][reduced-execution][sqe]") {
    constexpr auto chain_atom_count = std::size_t{12};
    for (const auto mode : {calculation::ExecutionSelectionKind::cutoff,
                            calculation::ExecutionSelectionKind::cover}) {
        CAPTURE(mode);

        const auto neutral = make_extended_components(0);
        const auto sqe = calculate_reduced(neutral, "sqe",
                                           {chargefw::test::make_sqe_ho_parameters("sqe")}, mode);
        check_component_total(sqe, 0, chain_atom_count, 0.0);
        check_component_total(sqe, chain_atom_count, chain_atom_count + 1, 0.0);

        const auto charged = make_extended_components(1);
        const auto sqeq0 = calculate_reduced(
            charged, "sqeq0", {chargefw::test::make_sqe_ho_parameters("sqeq0")}, mode);
        check_component_total(sqeq0, 0, chain_atom_count, 0.0);
        check_component_total(sqeq0, chain_atom_count, chain_atom_count + 1, 1.0);

        const auto sqeqp =
            calculate_reduced(charged, "sqeqp",
                              {chargefw::test::make_sqe_ho_parameters(
                                  "sqeqp", {.hydrogen_q0 = 0.25, .oxygen_q0 = -0.5})},
                              mode);
        check_component_total(sqeqp, 0, chain_atom_count, 0.0);
        check_component_total(sqeqp, chain_atom_count, chain_atom_count + 1, 1.0);
    }
}

TEST_CASE("SQE reduced execution fails when a fragment system is not positive definite",
          "[calculation][reduced-execution][sqe]") {
    const auto chain = make_extended_components(0);
    for (const auto mode : {calculation::ExecutionSelectionKind::cutoff,
                            calculation::ExecutionSelectionKind::cover}) {
        CAPTURE(mode);
        const auto result = calculate_application(calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::vector{chain}},
            .parameter_sets = {chargefw::test::make_sqe_ho_parameters("sqe", {.kappa = -50.0})},
            .method_id = "sqe",
            .execution_selection = calculation::ExecutionSelection{mode, 8.0}});

        CHECK(result.status == calculation::ExecutionStatus::numerical_failure);
        CHECK_FALSE(result.charges.has_value());
        REQUIRE(result.effective.has_value());
        CHECK(result.effective->method_id == "sqe");
        REQUIRE(result.failure_message.has_value());
        CHECK(result.failure_message->contains(mode == calculation::ExecutionSelectionKind::cutoff
                                                   ? "cutoff fragment around"
                                                   : "cover fragment around"));
    }
}

TEST_CASE("reduced approximation remains bounded across a truncated radius sweep",
          "[calculation][reduced-execution][accuracy]") {
    struct AccuracyCase {
        std::string_view method_id;
        core::Molecule molecule;
        std::vector<parameters::ParameterSet> parameter_sets;
        double max_rmsd;
        double max_mae;
        double max_maxabs;
    };
    const auto charged = make_extended_components(1);
    // Fixture-specific comparisons with full execution, not independent scientific reference
    // values or accuracy guarantees for other molecules and parameter sets.
    const auto cases = std::vector<AccuracyCase>{
        {.method_id = "abeem",
         .molecule = charged,
         .parameter_sets = {chargefw::test::make_abeem_ho_parameters()},
         .max_rmsd = 0.01,
         .max_mae = 0.01,
         .max_maxabs = 0.02},
        {.method_id = "eem",
         .molecule = charged,
         .parameter_sets = {chargefw::test::make_eem_parameters()},
         .max_rmsd = 0.03,
         .max_mae = 0.02,
         .max_maxabs = 0.08},
        {.method_id = "eqeq",
         .molecule = charged,
         .parameter_sets = {},
         .max_rmsd = 0.05,
         .max_mae = 0.03,
         .max_maxabs = 0.15},
        {.method_id = "eqeqc",
         .molecule = charged,
         .parameter_sets = {make_eqeqc_parameters()},
         .max_rmsd = 0.05,
         .max_mae = 0.03,
         .max_maxabs = 0.15},
        {.method_id = "qeq",
         .molecule = charged,
         .parameter_sets = {chargefw::test::make_qeq_ho_parameters()},
         .max_rmsd = 0.1,
         .max_mae = 0.05,
         .max_maxabs = 0.3},
        {.method_id = "sfkeem",
         .molecule = charged,
         .parameter_sets = {chargefw::test::make_sfkeem_ho_parameters()},
         .max_rmsd = 0.04,
         .max_mae = 0.02,
         .max_maxabs = 0.1},
        {.method_id = "sqe",
         .molecule = make_extended_components(0),
         .parameter_sets = {chargefw::test::make_sqe_ho_parameters("sqe")},
         .max_rmsd = 0.01,
         .max_mae = 0.01,
         .max_maxabs = 0.02},
        {.method_id = "sqeq0",
         .molecule = charged,
         .parameter_sets = {chargefw::test::make_sqe_ho_parameters("sqeq0")},
         .max_rmsd = 0.01,
         .max_mae = 0.01,
         .max_maxabs = 0.02},
        {.method_id = "sqeqp",
         .molecule = charged,
         .parameter_sets = {chargefw::test::make_sqe_ho_parameters(
             "sqeqp", {.hydrogen_q0 = 0.25, .oxygen_q0 = -0.5})},
         .max_rmsd = 0.01,
         .max_mae = 0.01,
         .max_maxabs = 0.02}};

    for (const auto& test_case : cases) {
        const features::PreparedMolecule prepared{test_case.molecule};
        const features::ConformerFeatures geometry{test_case.molecule, 0};
        const features::SpatialFragmentBuilder builder{prepared, geometry};
        CHECK(builder.build(0, 8.0).molecule().atom_count() < test_case.molecule.atom_count());

        const auto full =
            calculate_full(test_case.molecule, test_case.method_id, test_case.parameter_sets);
        for (const auto mode : {calculation::ExecutionSelectionKind::cutoff,
                                calculation::ExecutionSelectionKind::cover}) {
            for (const auto radius : {8.0, 12.0, 16.0}) {
                const auto reduced = calculate_reduced(test_case.molecule, test_case.method_id,
                                                       test_case.parameter_sets, mode, radius);
                const auto metrics = error_metrics(reduced, full);
                CAPTURE(test_case.method_id, mode, radius, metrics.rmsd, metrics.mae,
                        metrics.maxabs);
                CHECK(std::isfinite(metrics.rmsd));
                CHECK(std::isfinite(metrics.mae));
                CHECK(std::isfinite(metrics.maxabs));
                CHECK(metrics.rmsd < test_case.max_rmsd);
                CHECK(metrics.mae < test_case.max_mae);
                CHECK(metrics.maxabs < test_case.max_maxabs);
            }
        }
    }
}

TEST_CASE("SQE component totals coexist with intermolecular polarization",
          "[calculation][reduced-execution][sqe]") {
    const auto parameter_sets = std::vector{chargefw::test::make_sqe_ho_parameters("sqeq0")};
    const auto near = make_two_diatomic_components(4.0);
    const auto far = make_two_diatomic_components(30.0);

    const auto full_near = calculate_full(near, "sqeq0", parameter_sets);
    const auto full_far = calculate_full(far, "sqeq0", parameter_sets);
    check_component_total(full_near, 0, 2, 0.0);
    check_component_total(full_near, 2, 4, 0.0);
    CHECK(std::abs(full_near[0] - full_far[0]) > 1.0e-6);

    for (const auto mode : {calculation::ExecutionSelectionKind::cutoff,
                            calculation::ExecutionSelectionKind::cover}) {
        const auto reduced_near = calculate_reduced(near, "sqeq0", parameter_sets, mode);
        const auto reduced_far = calculate_reduced(far, "sqeq0", parameter_sets, mode);
        CAPTURE(mode);
        check_component_total(reduced_near, 0, 2, 0.0);
        check_component_total(reduced_near, 2, 4, 0.0);
        CHECK(std::abs(reduced_near[0] - reduced_far[0]) > 1.0e-6);
    }
}

TEST_CASE("reduced fragment calculations forward the complete fixed-source environment",
          "[calculation][reduced-execution][fixed-ions]") {
    const auto molecule = make_extended_components(1);
    const auto source = methods::FixedPointSource{core::Position{.x = 50.0, .y = 2.0}, 1.25};
    const auto sources = std::array{source};

    const auto exercise = [&](const std::string_view method_id,
                              const parameters::ParameterSet& parameter_set) {
        const features::PreparedMolecule prepared{molecule};
        const features::ConformerFeatures geometry{molecule, 0};
        const auto classification =
            parameters::classify_parameters(prepared.topology(), parameter_set);
        const auto* method = methods::method_registry().find(method_id);
        REQUIRE(method != nullptr);
        const auto options = methods::make_default_options(method->option_schema());
        const methods::ApplicableMethod selected{.method = method,
                                                 .parameter_set = &parameter_set,
                                                 .method_options = options,
                                                 .classifications = {classification}};
        const auto charge_context = calculation::detail::prepare_reduced_charge_context(
            selected, prepared, &classification);
        const features::SpatialFragmentBuilder builder{prepared, geometry};
        constexpr auto fragment_radius = 8.0;
        const auto fragment = builder.build(0, fragment_radius);
        REQUIRE(fragment.molecule().atom_count() < molecule.atom_count());
        CHECK(source.position.x > fragment_radius);

        const auto fragment_target =
            calculation::detail::fragment_target_charge(charge_context, fragment);

        const auto no_source_fragment = calculation::detail::calculate_fragment_charges(
            selected, &classification, fragment, charge_context);
        const auto zero_source = std::array{methods::FixedPointSource{source.position, 0.0}};
        const auto zero_source_fragment = calculation::detail::calculate_fragment_charges(
            selected, &classification, fragment, charge_context, zero_source);
        const auto sourced_fragment = calculation::detail::calculate_fragment_charges(
            selected, &classification, fragment, charge_context, sources);
        CHECK(std::abs(std::accumulate(no_source_fragment.values().begin(),
                                       no_source_fragment.values().end(), 0.0) -
                       fragment_target) < 1.0e-10);
        CHECK(std::abs(std::accumulate(sourced_fragment.values().begin(),
                                       sourced_fragment.values().end(), 0.0) -
                       fragment_target) < 1.0e-10);
        for (std::size_t atom_index = 0; atom_index < sourced_fragment.size(); ++atom_index) {
            CHECK(std::abs(zero_source_fragment[atom_index] - no_source_fragment[atom_index]) <
                  1.0e-14);
        }
        const auto center = fragment.center_local_atom_index();
        CHECK(std::abs(sourced_fragment[center] - no_source_fragment[center]) > 1.0e-9);

        const auto full_fragment = builder.build(0, 100.0);
        REQUIRE(full_fragment.molecule().atom_count() == molecule.atom_count());
        const auto no_source_whole = calculation::detail::calculate_fragment_charges(
            selected, &classification, full_fragment, charge_context);
        const auto sourced_whole = calculation::detail::calculate_fragment_charges(
            selected, &classification, full_fragment, charge_context, sources);
        const parameters::ParameterView parameter_view{parameter_set, classification};
        const auto full_target = core::total_formal_charge(molecule);
        const methods::CalculationInput full_input{prepared,  options,         full_target,
                                                   &geometry, &parameter_view, sources};
        const auto full_reference = method->calculate(full_input);
        const methods::CalculationInput no_source_input{prepared, options, full_target, &geometry,
                                                        &parameter_view};
        const auto no_source_reference = method->calculate(no_source_input);
        const auto source_indices = full_fragment.local_to_source_atom_indices();
        REQUIRE(source_indices.size() == molecule.atom_count());
        for (std::size_t local_index = 0; local_index < source_indices.size(); ++local_index) {
            const auto source_index = source_indices[local_index];
            CHECK(std::abs(sourced_whole[local_index] - full_reference[source_index]) < 1.0e-11);
            CHECK(std::abs(no_source_whole[local_index] - no_source_reference[source_index]) <
                  1.0e-11);
        }
    };

    exercise("eem", chargefw::test::make_eem_parameters());
    exercise("sqeqp", chargefw::test::make_sqe_ho_parameters(
                          "sqeqp", {.hydrogen_q0 = 0.25, .oxygen_q0 = -0.5}));
}

TEST_CASE("SQE+qp fragments start pieces of a split component from its reference charges",
          "[calculation][reduced-execution][sqeqp]") {
    // A U-shaped O-H chain whose arms run 2.5 Å apart and a separate H+: a 3 Å fragment around atom
    // 0 holds two unconnected pieces of the chain, atoms 0-2 and 18-19, and the H+.
    constexpr auto chain_atom_count = std::size_t{20};
    auto atoms = std::vector<core::Atom>{};
    auto bonds = std::vector<core::Bond>{};
    auto positions = std::vector<core::Position>{};
    for (std::size_t atom_index = 0; atom_index < chain_atom_count; ++atom_index) {
        atoms.emplace_back(atom_index % 2 == 0 ? 8 : 1, 0);
        const auto arm_index = atom_index < 10 ? atom_index : chain_atom_count - 1 - atom_index;
        positions.push_back(
            {.x = 1.4 * static_cast<double>(arm_index), .y = atom_index < 10 ? 0.0 : 2.5});
        if (atom_index != 0) {
            bonds.emplace_back(atom_index - 1, atom_index);
        }
    }
    atoms.emplace_back(1, 1);
    positions.push_back({.x = 1.0, .y = -2.0});
    const auto molecule = core::Molecule{
        std::move(atoms), std::move(bonds), {core::Conformer{std::move(positions)}}, "split-chain"};

    const auto parameter_set =
        chargefw::test::make_sqe_ho_parameters("sqeqp", {.hydrogen_q0 = 0.25, .oxygen_q0 = -0.5});
    const features::PreparedMolecule prepared{molecule};
    const features::ConformerFeatures geometry{molecule, 0};
    const auto classification = parameters::classify_parameters(prepared.topology(), parameter_set);
    const auto* method = methods::method_registry().find("sqeqp");
    REQUIRE(method != nullptr);
    const auto options = methods::make_default_options(method->option_schema());
    const methods::ApplicableMethod selected{.method = method,
                                             .parameter_set = &parameter_set,
                                             .method_options = options,
                                             .classifications = {classification}};
    const auto context =
        calculation::detail::prepare_reduced_charge_context(selected, prepared, &classification);
    const auto fragment = features::SpatialFragmentBuilder{prepared, geometry}.build(0, 3.0);
    const auto charges = calculation::detail::calculate_fragment_charges(selected, &classification,
                                                                         fragment, context);

    const auto source_indices = fragment.local_to_source_atom_indices();
    const auto pieces = features::connected_components(
        features::PreparedMolecule{fragment.molecule()}.topology().adjacency());
    REQUIRE(pieces.size() == 3);
    auto budget_differs_from_formal_charge = false;
    for (const auto& piece : pieces) {
        auto reference = 0.0;
        auto formal = 0.0;
        auto calculated = 0.0;
        for (const auto local_index : piece) {
            reference += context.reference_charges.at(source_indices[local_index]);
            formal += fragment.molecule().atom(local_index).formal_charge();
            calculated += charges[local_index];
        }
        CAPTURE(source_indices[piece.front()], reference, formal);
        CHECK(std::abs(calculated - reference) < 1.0e-12);
        budget_differs_from_formal_charge |= std::abs(reference - formal) > 0.1;
    }
    // Normalizing the pieces in the fragment would give them their formal charges instead.
    CHECK(budget_differs_from_formal_charge);
}

TEST_CASE("fixed-charge cutoff and cover diagnostics map active centers to original atoms",
          "[calculation][reduced-execution][fixed-ions]") {
    const auto original = core::MoleculeCollection{std::vector{
        core::Molecule{std::vector{core::Atom{12, 2}, core::Atom{1, 0}},
                       {},
                       {core::Conformer{{core::Position{.x = 0.0}, core::Position{.x = 20.0}}}},
                       "mapped-source"}}};
    const auto partition = calculation::detail::make_fixed_charge_partition(
        original, calculation::FixedIons{
                      .sources = {{.molecule_index = 0, .atom_index = 0, .charge = 0.25}}});
    const features::PreparedMoleculeCollection active{partition.active_molecules};
    const FailingFixedIonsFragmentMethod method;
    const methods::ApplicableMethod selected{.method = &method, .parameter_set = nullptr};

    for (const auto mode :
         {calculation::ExecutionMode::cutoff, calculation::ExecutionMode::cover}) {
        const auto policy = calculation::ExecutionPolicy{mode, 8.0};
        const auto calculate = [&] {
            if (mode == calculation::ExecutionMode::cutoff) {
                static_cast<void>(calculation::calculate_cutoff_charges(
                    selected, active, policy, 1, calculation::default_calculation_observer(),
                    &partition));
            } else {
                static_cast<void>(calculation::calculate_cover_charges(
                    selected, active, policy, 1, calculation::default_calculation_observer(),
                    &partition));
            }
        };
        const auto expected = mode == calculation::ExecutionMode::cutoff
                                  ? "cutoff fragment around source atom 2 failed"
                                  : "cover fragment around source atom 2 failed";
        CHECK_THROWS_MATCHES(calculate(), std::runtime_error,
                             snitch::matchers::with_what_contains{expected});
    }
}
