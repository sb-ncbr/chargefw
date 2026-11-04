#include <chargefw/core/molecule.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/methods/method_registry.h>

#include <snitch/snitch.hpp>

namespace methods = chargefw::methods;

TEST_CASE("VEEM rejects unsupported elements", "[methods][veem]") {
    const auto& registry = methods::method_registry();
    const auto* veem = registry.find("veem");

    REQUIRE(veem != nullptr);
    const chargefw::core::Molecule neon_molecule{{chargefw::core::Atom{10}}};
    const chargefw::features::PreparedMolecule prepared_neon{neon_molecule};
    const auto prerequisite_result = veem->check_method_prerequisites(
        {.prepared_molecule = prepared_neon, .method_options = {}});

    CHECK(!prerequisite_result);
    REQUIRE(prerequisite_result.issues().size() == 1);
    CHECK(prerequisite_result.issues()[0].kind ==
          methods::PrerequisiteIssueKind::unsupported_molecule);
    CHECK(prerequisite_result.issues()[0].location.atom_index == 0);
}
