#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/position.h>
#include <chargefw/features/conformer_features.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/features/spatial_fragment.h>

#include <algorithm>
#include <limits>
#include <ranges>
#include <stdexcept>
#include <utility>
#include <vector>

#include <snitch/snitch.hpp>

namespace core = chargefw::core;
namespace features = chargefw::features;

namespace {

auto make_molecule() -> core::Molecule {
    auto atoms = std::vector{core::Atom{6, 0, "C0"}, core::Atom{8, -1, "O1"},
                             core::Atom{7, 1, "N2"}, core::Atom{1, 0, "H3"}};
    auto bonds = std::vector{core::Bond{0, 1, core::BondOrder::DOUBLE},
                             core::Bond{1, 2, core::BondOrder::SINGLE},
                             core::Bond{2, 3, core::BondOrder::SINGLE}};
    auto conformers =
        std::vector{core::Conformer{{core::Position{.x = 0.0}, core::Position{.x = 1.0},
                                     core::Position{.x = 2.0}, core::Position{.x = 10.0}},
                                    "model-1"},
                    core::Conformer{{core::Position{.x = 0.0}, core::Position{.x = 2.0},
                                     core::Position{.x = 4.0}, core::Position{.x = 20.0}},
                                    "model-2"}};
    return core::Molecule{std::move(atoms), std::move(bonds), std::move(conformers), "fragment"};
}

} // namespace

TEST_CASE("spatial fragment builds from center atom and radius", "[features][spatial-fragment]") {
    const auto molecule = make_molecule();
    const features::PreparedMolecule prepared{molecule};
    const features::ConformerFeatures first_geometry{molecule, 0};
    const features::SpatialFragmentBuilder first_builder{prepared, first_geometry};
    const auto fragment = first_builder.build(1, 1.0);

    CHECK(fragment.molecule().name() == "fragment");
    CHECK(fragment.molecule().atom_count() == 3);
    CHECK(fragment.molecule().bond_count() == 2);
    CHECK(fragment.molecule().conformer_count() == 1);
    CHECK(fragment.molecule().conformer(0).name() == "model-1");
    CHECK(fragment.source_conformer_index() == 0);
    CHECK(fragment.center_local_atom_index() == 1);

    CHECK(std::ranges::equal(fragment.local_to_source_atom_indices(),
                             std::vector<std::size_t>{0, 1, 2}));
    CHECK(std::ranges::equal(fragment.local_to_source_bond_indices(),
                             std::vector<std::size_t>{0, 1}));

    CHECK(fragment.molecule().atom(1).atomic_number() == 8);
    CHECK(fragment.molecule().atom(1).formal_charge() == -1);
    CHECK(fragment.molecule().atom(1).name() == "O1");
    CHECK(fragment.molecule().bond(0).first_atom_index() == 0);
    CHECK(fragment.molecule().bond(0).second_atom_index() == 1);
    CHECK(fragment.molecule().bond(0).order() == core::BondOrder::DOUBLE);
    CHECK(fragment.molecule().bond(1).first_atom_index() == 1);
    CHECK(fragment.molecule().bond(1).second_atom_index() == 2);
    CHECK(fragment.molecule().conformer(0)[0].x == 0.0);
    CHECK(fragment.molecule().conformer(0)[1].x == 1.0);
    CHECK(fragment.molecule().conformer(0)[2].x == 2.0);
}

TEST_CASE("spatial fragment maps a radius selection to local indices",
          "[features][spatial-fragment]") {
    const auto molecule = make_molecule();
    const features::PreparedMolecule prepared{molecule};
    const features::ConformerFeatures first_geometry{molecule, 0};
    const features::SpatialFragmentBuilder first_builder{prepared, first_geometry};
    const auto fragment = first_builder.build(2, 1.0);

    REQUIRE(fragment.molecule().atom_count() == 2);
    REQUIRE(fragment.molecule().bond_count() == 1);
    CHECK(std::ranges::equal(fragment.local_to_source_atom_indices(),
                             std::vector<std::size_t>{1, 2}));
    CHECK(std::ranges::equal(fragment.local_to_source_bond_indices(), std::vector<std::size_t>{1}));
    CHECK(fragment.center_local_atom_index() == 1);
    CHECK(fragment.molecule().bond(0).first_atom_index() == 0);
    CHECK(fragment.molecule().bond(0).second_atom_index() == 1);
    CHECK(fragment.molecule().atom(0).name() == "O1");
    CHECK(fragment.molecule().atom(1).name() == "N2");
    CHECK(fragment.molecule().conformer(0)[0].x == 1.0);
    CHECK(fragment.molecule().conformer(0)[1].x == 2.0);
}

TEST_CASE("spatial fragment includes all atoms within boundary radius",
          "[features][spatial-fragment]") {
    const auto molecule = make_molecule();
    const features::PreparedMolecule prepared{molecule};
    const features::ConformerFeatures first_geometry{molecule, 0};
    const features::SpatialFragmentBuilder first_builder{prepared, first_geometry};
    const auto boundary_fragment = first_builder.build(2, 8.0);

    CHECK(boundary_fragment.molecule().atom_count() == 4);
    CHECK(boundary_fragment.center_local_atom_index() == 2);
    CHECK(std::ranges::equal(boundary_fragment.local_to_source_atom_indices(),
                             std::vector<std::size_t>{0, 1, 2, 3}));
}

TEST_CASE("spatial fragment uses second conformer geometry", "[features][spatial-fragment]") {
    const auto molecule = make_molecule();
    const features::PreparedMolecule prepared{molecule};
    const features::ConformerFeatures second_geometry{molecule, 1};
    const features::SpatialFragmentBuilder second_builder{prepared, second_geometry};
    const auto second_conformer_fragment = second_builder.build(1, 2.0);

    CHECK(second_conformer_fragment.molecule().conformer(0).name() == "model-2");
    CHECK(second_conformer_fragment.molecule().conformer(0)[2].x == 4.0);
}

TEST_CASE("spatial fragment builder is deterministic for identical parameters",
          "[features][spatial-fragment]") {
    const auto molecule = make_molecule();
    const features::PreparedMolecule prepared{molecule};
    const features::ConformerFeatures first_geometry{molecule, 0};
    const features::SpatialFragmentBuilder first_builder{prepared, first_geometry};
    const auto fragment = first_builder.build(1, 1.0);
    const auto repeated_fragment = first_builder.build(1, 1.0);

    CHECK(std::ranges::equal(repeated_fragment.local_to_source_atom_indices(),
                             fragment.local_to_source_atom_indices()));
}

TEST_CASE("spatial fragment supports single-atom molecules", "[features][spatial-fragment]") {
    const core::Molecule single_atom_molecule{
        std::vector{core::Atom{1, 0, "H"}},
        {},
        std::vector{core::Conformer{{core::Position{.x = 3.0}}, "single"}},
        "single"};
    const features::PreparedMolecule single_atom_prepared{single_atom_molecule};
    const features::ConformerFeatures single_atom_geometry{single_atom_molecule, 0};
    const features::SpatialFragmentBuilder single_atom_builder{single_atom_prepared,
                                                               single_atom_geometry};
    const auto single_atom_fragment = single_atom_builder.build(0, 1.0);

    CHECK(single_atom_fragment.molecule().atom_count() == 1);
    CHECK(single_atom_fragment.molecule().bond_count() == 0);
    CHECK(single_atom_fragment.center_local_atom_index() == 0);
}

TEST_CASE("spatial fragment builder rejects mismatched geometry and prepared molecule",
          "[features][spatial-fragment]") {
    const auto molecule = make_molecule();
    const features::PreparedMolecule prepared{molecule};

    const core::Molecule other_molecule{
        std::vector{core::Atom{1, 0, "H"}},
        {},
        std::vector{core::Conformer{{core::Position{.x = 3.0}}, "single"}},
        "single"};
    const features::ConformerFeatures unrelated_geometry{other_molecule, 0};

    CHECK_THROWS_AS((features::SpatialFragmentBuilder{prepared, unrelated_geometry}),
                    std::invalid_argument);
}

TEST_CASE("spatial fragment builder rejects non-finite conformer geometry",
          "[features][spatial-fragment]") {
    const core::Molecule nonfinite_molecule{
        std::vector{core::Atom{1, 0, "H0"}, core::Atom{1, 0, "H1"}},
        {},
        std::vector{core::Conformer{
            {core::Position{}, core::Position{.x = std::numeric_limits<double>::quiet_NaN()}},
            "nonfinite"}},
        "nonfinite"};
    const features::PreparedMolecule nonfinite_prepared{nonfinite_molecule};
    const features::ConformerFeatures nonfinite_geometry{nonfinite_molecule, 0};

    CHECK_THROWS_AS((features::SpatialFragmentBuilder{nonfinite_prepared, nonfinite_geometry}),
                    std::invalid_argument);
}

TEST_CASE("spatial fragment builder rejects invalid center, radius, and empty molecule",
          "[features][spatial-fragment]") {
    const auto molecule = make_molecule();
    const features::PreparedMolecule prepared{molecule};
    const features::ConformerFeatures first_geometry{molecule, 0};
    const features::SpatialFragmentBuilder first_builder{prepared, first_geometry};

    const core::Molecule empty_molecule{{}, {}, {core::Conformer{{}, "empty"}}, "empty"};
    const features::PreparedMolecule empty_prepared{empty_molecule};
    const features::ConformerFeatures empty_geometry{empty_molecule, 0};
    const features::SpatialFragmentBuilder empty_builder{empty_prepared, empty_geometry};

    CHECK_THROWS_AS(empty_builder.build(0, 1.0), std::out_of_range);
    CHECK_THROWS_AS(first_builder.build(4, 1.0), std::out_of_range);
    CHECK_THROWS_AS(first_builder.build(1, 0.0), std::invalid_argument);
    CHECK_THROWS_AS(first_builder.build(1, -1.0), std::invalid_argument);
    CHECK_THROWS_AS(first_builder.build(1, std::numeric_limits<double>::quiet_NaN()),
                    std::invalid_argument);
    CHECK_THROWS_AS(first_builder.build(1, std::numeric_limits<double>::infinity()),
                    std::invalid_argument);
}
