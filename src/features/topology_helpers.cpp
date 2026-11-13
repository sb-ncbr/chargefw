#include "features/topology_helpers.h"

#include <chargefw/features/topology_features.h>

#include <algorithm>
#include <queue>
#include <vector>

namespace chargefw::features {
namespace {

auto breadth_first_distances(const TopologyFeatures& topology, const std::size_t source)
    -> std::vector<int> {
    auto distances = std::vector(topology.molecule().atom_count(), -1);
    auto queue = std::queue<std::size_t>{};

    distances[source] = 0;
    queue.push(source);

    while (!queue.empty()) {
        const auto current = queue.front();
        queue.pop();

        const auto current_distance = distances[current];

        for (const auto neighbor : topology.neighbor_indices(current)) {
            if (distances[neighbor] != -1) {
                continue;
            }

            distances[neighbor] = current_distance + 1;
            queue.push(neighbor);
        }
    }

    return distances;
}

} // namespace

auto all_pairs_bond_distances(const TopologyFeatures& topology) -> std::vector<std::vector<int>> {
    const auto atom_count = topology.molecule().atom_count();
    auto distances = std::vector<std::vector<int>>{};
    distances.reserve(atom_count);

    for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
        distances.push_back(breadth_first_distances(topology, atom_index));
    }

    return distances;
}

auto connected_components(const TopologyFeatures& topology)
    -> std::vector<std::vector<std::size_t>> {
    const auto atom_count = topology.molecule().atom_count();
    auto components = std::vector<std::vector<std::size_t>>{};
    auto visited = std::vector<bool>(atom_count, false);
    auto pending = std::vector<std::size_t>{};

    for (std::size_t first_atom_index = 0; first_atom_index < atom_count; ++first_atom_index) {
        if (visited[first_atom_index]) {
            continue;
        }

        auto& component = components.emplace_back();
        visited[first_atom_index] = true;
        pending.push_back(first_atom_index);

        while (!pending.empty()) {
            const auto atom_index = pending.back();
            pending.pop_back();
            component.push_back(atom_index);

            for (const auto neighbor_index : topology.neighbor_indices(atom_index)) {
                if (!visited[neighbor_index]) {
                    visited[neighbor_index] = true;
                    pending.push_back(neighbor_index);
                }
            }
        }
    }

    return components;
}

auto is_connected(const TopologyFeatures& topology) -> bool {
    if (topology.molecule().atom_count() == 0) {
        return true;
    }

    const auto distances = breadth_first_distances(topology, 0);

    return std::ranges::all_of(distances, [](const int distance) -> bool { return distance >= 0; });
}

} // namespace chargefw::features
