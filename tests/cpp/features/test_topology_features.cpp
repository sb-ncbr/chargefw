#include "support/test_molecules.h"

#include "features/topology_helpers.h"

#include <chargefw/core/molecule.h>
#include <chargefw/features/topology_features.h>

#include <span>
#include <stdexcept>
#include <vector>

#include <snitch/snitch.hpp>

namespace core = chargefw::core;
namespace features = chargefw::features;

namespace {

auto to_vector(std::span<const std::size_t> values) -> std::vector<std::size_t> {
    return {values.begin(), values.end()};
}

} // namespace

TEST_CASE("topology features expose molecule reference and graph degree", "[features][topology]") {
    const auto water = chargefw::test::make_water();
    const features::TopologyFeatures topology{water};

    CHECK(&topology.molecule() == &water);
    CHECK(topology.degree(0) == 2);
    CHECK(topology.degree(1) == 1);
    CHECK(topology.degree(2) == 1);
}

TEST_CASE("topology features expose neighbor indices", "[features][topology]") {
    const auto water = chargefw::test::make_water();
    const features::TopologyFeatures topology{water};

    const auto oxygen_neighbors = to_vector(topology.neighbor_indices(0));
    CHECK(oxygen_neighbors == std::vector<std::size_t>{1, 2});

    const auto h1_neighbors = to_vector(topology.neighbor_indices(1));
    CHECK(h1_neighbors == std::vector<std::size_t>{0});
}

TEST_CASE("topology features expose incident bond indices", "[features][topology]") {
    const auto water = chargefw::test::make_water();
    const features::TopologyFeatures topology{water};

    CHECK(to_vector(topology.incident_bond_indices(0)) == std::vector<std::size_t>{0, 1});
    CHECK(to_vector(topology.incident_bond_indices(1)) == std::vector<std::size_t>{0});
}

TEST_CASE("topology features reject invalid atom indices", "[features][topology]") {
    const auto water = chargefw::test::make_water();
    const features::TopologyFeatures topology{water};

    CHECK_THROWS_AS(topology.degree(3), std::out_of_range);
}

TEST_CASE("topology helpers identify connected components", "[features][topology]") {
    const core::Molecule pieces{
        {core::Atom{1}, core::Atom{1}, core::Atom{1}, core::Atom{1}, core::Atom{2}},
        {core::Bond{0, 1}, core::Bond{2, 3}}};
    const core::Molecule empty{{}};

    const auto components = features::connected_components(features::TopologyFeatures{pieces});

    REQUIRE(components.size() == 3);
    CHECK(components[0] == std::vector<std::size_t>{0, 1});
    CHECK(components[1] == std::vector<std::size_t>{2, 3});
    CHECK(components[2] == std::vector<std::size_t>{4});
    CHECK(features::connected_components(features::TopologyFeatures{empty}).empty());
}
