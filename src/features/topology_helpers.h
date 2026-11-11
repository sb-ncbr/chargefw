#pragma once

#include <cstddef>
#include <vector>

namespace chargefw::features {

auto all_pairs_bond_distances(const std::vector<std::vector<std::size_t>>& adjacency)
    -> std::vector<std::vector<int>>;

auto connected_components(const std::vector<std::vector<std::size_t>>& adjacency)
    -> std::vector<std::vector<std::size_t>>;

auto is_connected(const std::vector<std::vector<std::size_t>>& adjacency) -> bool;

} // namespace chargefw::features
