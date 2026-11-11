#include "calculation/execution_assessment.h"
#include "support/test_methods.h"
#include "support/test_molecules.h"

#include <chargefw/calculation/execution_issues.h>
#include <chargefw/calculation/execution_policy.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method_registry.h>
#include <chargefw/methods/method_requirements.h>

#include <cstddef>
#include <optional>
#include <snitch/snitch.hpp>
#include <vector>

namespace calculation = chargefw::calculation;
namespace core = chargefw::core;
namespace features = chargefw::features;
namespace methods = chargefw::methods;

using calculation::ExecutionMode;
using calculation::detail::assess_execution;
using calculation::detail::ExecutionAvailability;

TEST_CASE("full execution warns above the cutoff threshold only for expensive methods",
          "[calculation][execution-assessment]") {
    // Three and two atoms.
    const core::MoleculeCollection collection{
        std::vector{chargefw::test::make_water(), chargefw::test::make_formally_charged_pair()}};
    const features::PreparedMoleculeCollection molecules{collection};

    const chargefw::test::StubMethod inexpensive{"resource-test"};
    CHECK(
        assess_execution(inexpensive, molecules, {.cutoff_atom_threshold = 1}, ExecutionMode::full)
            .availability == ExecutionAvailability::available);

    const chargefw::test::StubMethod expensive{
        "resource-test",
        {.resources = {.time = methods::ComplexityTerm::atoms_cubed,
                       .memory = methods::ComplexityTerm::atoms_squared}}};
    const auto full = [&](const calculation::ResourcePolicy& policy) {
        return assess_execution(expensive, molecules, policy, ExecutionMode::full);
    };

    const auto below = full({.cutoff_atom_threshold = 3});
    CHECK(below.availability == ExecutionAvailability::available);
    CHECK(below.issues.empty());

    const auto above = full({.cutoff_atom_threshold = 2});
    CHECK(above.availability == ExecutionAvailability::available_with_warning);
    REQUIRE(above.issues.size() == 1);
    CHECK(above.issues[0].kind == calculation::ExecutionIssueKind::resource_threshold_exceeded);
    CHECK(above.issues[0].molecule_index == std::optional<std::size_t>{0});

    // One issue per molecule above the threshold; none without a threshold.
    CHECK(full({.cutoff_atom_threshold = 1}).issues.size() == molecules.size());
    const auto unlimited =
        full({.cutoff_atom_threshold = std::nullopt, .cover_atom_threshold = std::nullopt});
    CHECK(unlimited.availability == ExecutionAvailability::available);
    CHECK(unlimited.issues.empty());

    // A bundled method with expensive full execution and no reduced execution.
    const auto* mgc = methods::method_registry().find("mgc");
    REQUIRE(mgc != nullptr);
    const core::MoleculeCollection graph{std::vector{chargefw::test::make_water_graph()}};
    const features::PreparedMoleculeCollection prepared_graph{graph};
    const calculation::ResourcePolicy policy{.cutoff_atom_threshold = 2};
    CHECK(assess_execution(*mgc, prepared_graph, policy, ExecutionMode::full).availability ==
          ExecutionAvailability::available_with_warning);
    for (const auto mode : {ExecutionMode::cutoff, ExecutionMode::cover}) {
        const auto reduced = assess_execution(*mgc, prepared_graph, policy, mode);
        CHECK(reduced.availability == ExecutionAvailability::unsupported);
        REQUIRE(reduced.issues.size() == 1);
        CHECK(reduced.issues[0].kind ==
              calculation::ExecutionIssueKind::unsupported_execution_mode);
    }
}

TEST_CASE("reduced execution requires coordinates and declared support",
          "[calculation][execution-assessment]") {
    const core::MoleculeCollection collection{std::vector{chargefw::test::make_water()}};
    const features::PreparedMoleculeCollection molecules{collection};

    const chargefw::test::StubMethod topology_only{
        "resource-test", {.resources = {.supports_cutoff = true, .supports_cover = true}}};
    const chargefw::test::StubMethod spatial{
        "resource-test",
        {.coordinates = true,
         .resources = {.supports_cutoff = true,
                       .supports_cover = true,
                       .reduced_charge_policy =
                           methods::ReducedChargePolicy::uniform_target_global}}};
    for (const auto mode : {ExecutionMode::cutoff, ExecutionMode::cover}) {
        CHECK(assess_execution(topology_only, molecules, {}, mode).availability ==
              ExecutionAvailability::unsupported);
        CHECK(assess_execution(spatial, molecules, {}, mode).availability ==
              ExecutionAvailability::available);
    }

    // Cutoff execution of a molecule above the cover threshold warns; cover does not.
    const calculation::ResourcePolicy small_cover{.cutoff_atom_threshold = 1,
                                                  .cover_atom_threshold = 2};
    CHECK(assess_execution(spatial, molecules, small_cover, ExecutionMode::cutoff).availability ==
          ExecutionAvailability::available_with_warning);
    CHECK(assess_execution(spatial, molecules, small_cover, ExecutionMode::cover).availability ==
          ExecutionAvailability::available);
}
