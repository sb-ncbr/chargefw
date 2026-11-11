#include "features/topology_helpers.h"

#include <algorithm>
#include <queue>
#include <vector>

namespace chargefw::features {
namespace {

auto breadth_first_distances(const std::vector<std::vector<std::size_t>>& adjacency,
                             const std::size_t source) -> std::vector<int> {
    auto distances = std::vector(adjacency.size(), -1);
    auto queue = std::queue<std::size_t>{};

    distances[source] = 0;
    queue.push(source);

    while (!queue.empty()) {
        const auto current = queue.front();
        queue.pop();

        const auto current_distance = distances[current];

        for (const auto neighbor : adjacency[current]) {
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

auto all_pairs_bond_distances(const std::vector<std::vector<std::size_t>>& adjacency)
    -> std::vector<std::vector<int>> {
    auto distances = std::vector<std::vector<int>>{};
    distances.reserve(adjacency.size());

    for (std::size_t atom_index = 0; atom_index < adjacency.size(); ++atom_index) {
        distances.push_back(breadth_first_distances(adjacency, atom_index));
    }

    return distances;
}

auto connected_components(const std::vector<std::vector<std::size_t>>& adjacency)
    -> std::vector<std::vector<std::size_t>> {
    auto components = std::vector<std::vector<std::size_t>>{};
    auto visited = std::vector<bool>(adjacency.size(), false);
    auto pending = std::vector<std::size_t>{};

    for (std::size_t first_atom_index = 0; first_atom_index < adjacency.size();
         ++first_atom_index) {
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

            for (const auto neighbor_index : adjacency[atom_index]) {
                if (!visited[neighbor_index]) {
                    visited[neighbor_index] = true;
                    pending.push_back(neighbor_index);
                }
            }
        }
    }

    return components;
}

auto is_connected(const std::vector<std::vector<std::size_t>>& adjacency) -> bool {
    if (adjacency.empty()) {
        return true;
    }

    const auto distances = breadth_first_distances(adjacency, 0);

    return std::ranges::all_of(distances, [](const int distance) -> bool { return distance >= 0; });
}

} // namespace chargefw::features
