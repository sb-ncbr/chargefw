#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/molecule.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/methods/method_registry.h>

#include <snitch/snitch.hpp>

namespace methods = chargefw::methods;

TEST_CASE("TSEF rejects disconnected molecules", "[methods][tsef]") {
    const auto* tsef = methods::method_registry().find("tsef");
    REQUIRE(tsef != nullptr);

    // Bonded, so the general bond requirement passes, but in two parts.
    const chargefw::core::Molecule disconnected{
        {chargefw::core::Atom{1}, chargefw::core::Atom{8}, chargefw::core::Atom{1}},
        {chargefw::core::Bond{0, 1, chargefw::core::BondOrder::SINGLE}}};
    const chargefw::features::PreparedMolecule prepared_pair{disconnected};
    const methods::MethodOptions options;
    const auto prerequisite_result = tsef->check_method_prerequisites(
        {.prepared_molecule = prepared_pair, .method_options = options});

    CHECK(!prerequisite_result);
    REQUIRE(prerequisite_result.issues().size() == 1);
    CHECK(prerequisite_result.issues()[0].kind ==
          methods::PrerequisiteIssueKind::unsupported_molecule);
}
