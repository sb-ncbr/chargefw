#include "support/test_molecules.h"

#include <chargefw/core/atom.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/core/position.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/features/prepared_molecule_collection.h>

#include <cstddef>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

#include <snitch/snitch.hpp>

namespace core = chargefw::core;
namespace features = chargefw::features;

static_assert(!std::is_constructible_v<features::PreparedMolecule, core::Molecule&&>);
static_assert(std::is_copy_constructible_v<features::PreparedMolecule> &&
              std::is_copy_assignable_v<features::PreparedMolecule> &&
              std::is_nothrow_move_constructible_v<features::PreparedMolecule> &&
              std::is_nothrow_move_assignable_v<features::PreparedMolecule>);
static_assert(std::is_copy_constructible_v<features::PreparedMoleculeCollection>);
static_assert(
    !std::is_constructible_v<features::PreparedMoleculeCollection, core::MoleculeCollection&&>);

TEST_CASE("prepared molecules retain source identity and cached topology", "[features][prepared]") {
    const auto water = chargefw::test::make_water_graph();
    const features::PreparedMolecule prepared{water};

    CHECK(&prepared.molecule() == &water);
    CHECK(&prepared.topology().molecule() == &water);
    CHECK(prepared.topology().degree(0) == 2);
}

TEST_CASE("prepared molecules find geometry defects once per conformer", "[features][prepared]") {
    const auto infinity = std::numeric_limits<double>::infinity();
    // Carbon, hydrogen, and chlorine have covalent radii of 0.76, 0.31, and 1.02 Å.
    const core::Molecule molecule{
        std::vector{core::Atom{6}, core::Atom{1}, core::Atom{17}},
        {},
        std::vector{
            core::Conformer{{core::Position{}, core::Position{.x = 1.0}, core::Position{.x = 1.0}},
                            "coincident"},
            core::Conformer{
                {core::Position{}, core::Position{.x = infinity}, core::Position{.x = 2.0}},
                "nonfinite"},
            // Identical positions are not looked for when a coordinate is not finite.
            core::Conformer{{core::Position{}, core::Position{}, core::Position{.x = infinity}},
                            "both"},
            // The hydrogen is the nearest neighbor of the carbon, but only the chlorine is closer
            // than 0.7 times the sum of the radii: 1.0 / 1.78 = 0.56, against 0.9 / 1.07 = 0.84.
            core::Conformer{{core::Position{}, core::Position{.x = 0.9}, core::Position{.y = 1.0}},
                            "close"},
            // Close atoms are not reported when two atoms are at identical positions; -0.0 and
            // 0.0 are the same position.
            core::Conformer{{core::Position{}, core::Position{.x = -0.0}, core::Position{.y = 0.5}},
                            "coincident and close"}}};
    const features::PreparedMolecule prepared{molecule};

    const auto defects = prepared.geometry_defects();
    REQUIRE(defects.size() == 5);
    CHECK_FALSE(defects[0].nonfinite_atom.has_value());
    CHECK(defects[0].coincident_atoms == std::pair<std::size_t, std::size_t>{1, 2});
    CHECK_FALSE(defects[0].close_atoms.has_value());
    CHECK(defects[1].nonfinite_atom == 1);
    CHECK_FALSE(defects[1].coincident_atoms.has_value());
    CHECK(defects[2].nonfinite_atom == 2);
    CHECK_FALSE(defects[2].coincident_atoms.has_value());
    CHECK_FALSE(defects[3].coincident_atoms.has_value());
    REQUIRE(defects[3].close_atoms.has_value());
    CHECK(defects[3].close_atoms->closest == std::pair<std::size_t, std::size_t>{0, 2});
    CHECK(defects[3].close_atoms->distance == 1.0);
    CHECK(defects[3].close_atoms->pair_count == 1);
    CHECK(defects[4].coincident_atoms == std::pair<std::size_t, std::size_t>{0, 1});
    CHECK_FALSE(defects[4].close_atoms.has_value());
    CHECK(prepared.geometry_defects().data() == defects.data());

    // Copies share the defects, and geometry queries remain valid on a moved-from object.
    // NOLINTNEXTLINE(performance-unnecessary-copy-initialization): the copy is under test.
    const auto copy = prepared;
    CHECK(&copy.molecule() == &molecule);
    CHECK(copy.geometry_defects().data() == defects.data());
    auto source = features::PreparedMolecule{molecule};
    const auto moved = std::move(source);
    CHECK(moved.geometry_defects().size() == 5);
    // NOLINTNEXTLINE(bugprone-use-after-move): the moved-from source is under test.
    CHECK(source.geometry_defects().size() == 5);
}

TEST_CASE("prepared molecule collection retains each source molecule", "[features][prepared]") {
    const auto collection = core::MoleculeCollection{
        std::vector{chargefw::test::make_water_graph(), chargefw::test::make_water_graph()},
        "waters"};
    const features::PreparedMoleculeCollection prepared{collection};

    REQUIRE(prepared.size() == collection.size());
    CHECK(&prepared[0].molecule() == &collection[0]);
    CHECK(&prepared[1].molecule() == &collection[1]);
    CHECK(prepared[0].topology().degree(0) == 2);
    CHECK(prepared[1].topology().degree(0) == 2);
}
