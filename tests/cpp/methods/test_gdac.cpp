#include "support/test_molecules.h"

#include <chargefw/core/atom.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/methods/method_registry.h>

#include <snitch/snitch.hpp>

namespace features = chargefw::features;
namespace methods = chargefw::methods;

TEST_CASE("GDAC rejects missing features and unsupported elements", "[methods][gdac]") {
    const auto& registry = methods::method_registry();
    const auto* gdac = registry.find("gdac");

    REQUIRE(gdac != nullptr);
    const auto options = methods::make_default_options(gdac->option_schema());

    const auto charged_pair = chargefw::test::make_formally_charged_pair();
    const features::PreparedMolecule prepared_charged_pair{charged_pair};

    const auto prerequisite_result = gdac->check_method_prerequisites(
        {.prepared_molecule = prepared_charged_pair, .method_options = options});

    CHECK(!prerequisite_result);
    REQUIRE(!prerequisite_result.issues.empty());
    CHECK(prerequisite_result.issues[0].kind == methods::PrerequisiteIssueKind::missing_feature);

    const chargefw::core::Molecule rubidium_molecule{
        {chargefw::core::Atom{37}}, {}, {chargefw::core::Conformer{{chargefw::core::Position{}}}}};
    const features::PreparedMolecule prepared_rubidium{rubidium_molecule};
    const auto rubidium_prerequisite_result = gdac->check_method_prerequisites(
        {.prepared_molecule = prepared_rubidium, .method_options = options});

    CHECK(!rubidium_prerequisite_result);
    REQUIRE(rubidium_prerequisite_result.issues.size() == 1);
    CHECK(rubidium_prerequisite_result.issues[0].kind ==
          methods::PrerequisiteIssueKind::unsupported_molecule);
    CHECK(rubidium_prerequisite_result.issues[0].location.atom_index == 0);
}
