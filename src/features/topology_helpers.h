#pragma once

#include <cstddef>
#include <vector>

namespace chargefw::features {

class TopologyFeatures;

auto all_pairs_bond_distances(const TopologyFeatures& topology) -> std::vector<std::vector<int>>;

auto connected_components(const TopologyFeatures& topology)
    -> std::vector<std::vector<std::size_t>>;

auto is_connected(const TopologyFeatures& topology) -> bool;

} // namespace chargefw::features
