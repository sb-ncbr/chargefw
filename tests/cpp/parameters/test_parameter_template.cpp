#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/parameters/classification/parameter_classifier.h>
#include <chargefw/parameters/classification/parameter_template.h>
#include <chargefw/parameters/models/parameter_key.h>
#include <chargefw/parameters/models/parameter_values.h>

#include <array>
#include <cstddef>
#include <snitch/snitch.hpp>
#include <utility>
#include <vector>

namespace core = chargefw::core;
namespace features = chargefw::features;
namespace parameters = chargefw::parameters;
using parameters::AtomParameterClassificationKind;
using parameters::BondParameterClassificationKind;

namespace {

auto make_formaldehyde() -> core::Molecule {
    return core::Molecule{{core::Atom{6, 0, "C"}, core::Atom{8, 0, "O"}, core::Atom{1, 0, "H1"},
                           core::Atom{1, 0, "H2"}},
                          {core::Bond{0, 1, core::BondOrder::DOUBLE},
                           core::Bond{0, 2, core::BondOrder::SINGLE},
                           core::Bond{3, 0, core::BondOrder::SINGLE}},
                          {},
                          "formaldehyde"};
}

auto make_collection() -> core::MoleculeCollection {
    return core::MoleculeCollection{
        std::vector{chargefw::test::make_water_graph(), make_formaldehyde()}, "training"};
}

const auto sqe_names =
    parameters::ParameterNames{.atom = {"electronegativity", "hardness"}, .bond = {"kappa"}};

} // namespace

TEST_CASE("parameter template has one zero entry per occurring type, in key order",
          "[parameters][parameter-template]") {
    const auto collection = make_collection();
    const auto prepared = features::PreparedMoleculeCollection{collection};

    const auto set = parameters::make_parameter_template(
        prepared, sqe_names, AtomParameterClassificationKind::HIGHEST_BOND_ORDER,
        BondParameterClassificationKind::BOND_ORDER, {.id = "template", .method_id = "sqe"});

    using chargefw::test::hbo_atom_key;
    const auto expected_atoms = std::vector{hbo_atom_key(1, "1"), hbo_atom_key(6, "2"),
                                            hbo_atom_key(8, "1"), hbo_atom_key(8, "2")};
    REQUIRE(set.atom().size() == expected_atoms.size());
    for (std::size_t index = 0; index < expected_atoms.size(); ++index) {
        CHECK(set.atom()[index].key == expected_atoms[index]);
    }

    const auto bond_key = [](parameters::AtomParameterKey first,
                             parameters::AtomParameterKey second, std::string type) {
        return parameters::BondParameterKey{
            .first_atom = std::move(first),
            .second_atom = std::move(second),
            .bond = {.classification = BondParameterClassificationKind::BOND_ORDER,
                     .type = std::move(type)}};
    };
    // Each bond appears once, with the smaller atom key first, whatever its order in the molecule.
    const auto expected_bonds =
        std::vector{bond_key(hbo_atom_key(1, "1"), hbo_atom_key(6, "2"), "1"),
                    bond_key(hbo_atom_key(1, "1"), hbo_atom_key(8, "1"), "1"),
                    bond_key(hbo_atom_key(6, "2"), hbo_atom_key(8, "2"), "2")};
    REQUIRE(set.bond().size() == expected_bonds.size());
    for (std::size_t index = 0; index < expected_bonds.size(); ++index) {
        CHECK(set.bond()[index].key == expected_bonds[index]);
    }

    CHECK(set.method_id() == "sqe");
    CHECK(set.common().empty());
    CHECK(parameters::parameter_values(set) == std::vector<double>(4 * 2 + 3 * 1, 0.0));
}

TEST_CASE("parameter template matches every atom and bond of its molecules",
          "[parameters][parameter-template]") {
    const auto collection = make_collection();
    const auto prepared = features::PreparedMoleculeCollection{collection};

    for (const auto& [atom_types, bond_types] :
         std::array{std::pair{AtomParameterClassificationKind::PLAIN,
                              BondParameterClassificationKind::PLAIN},
                    std::pair{AtomParameterClassificationKind::HIGHEST_BOND_ORDER,
                              BondParameterClassificationKind::BOND_ORDER},
                    std::pair{AtomParameterClassificationKind::BONDED_ELEMENTS,
                              BondParameterClassificationKind::BOND_ORDER}}) {
        CAPTURE(parameters::to_string(atom_types), parameters::to_string(bond_types));
        const auto set = parameters::make_parameter_template(
            prepared, sqe_names, atom_types, bond_types, {.id = "template", .method_id = "sqe"});

        for (const auto& molecule : prepared.molecules()) {
            CAPTURE(molecule.molecule().name());
            const auto result =
                parameters::try_classify_parameters(molecule.molecule(), molecule.topology(), set);
            CHECK(static_cast<bool>(result));
        }
    }
}

TEST_CASE("parameter template makes only the sections it has names for",
          "[parameters][parameter-template]") {
    const auto collection = make_collection();
    const auto prepared = features::PreparedMoleculeCollection{collection};

    const auto set = parameters::make_parameter_template(
        prepared, {.common = {"kappa"}, .atom = {"A", "B"}}, AtomParameterClassificationKind::PLAIN,
        BondParameterClassificationKind::PLAIN, {.id = "template", .method_id = "eem"});

    CHECK(set.common().size() == 1);
    CHECK(set.atom().size() == 3); // H, C, O
    CHECK(set.bond().empty());
}
