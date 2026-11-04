#include <chargefw/core/molecule.h>
#include <chargefw/core/periodic_table.h>
#include <chargefw/features/conformer_features.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/features/topology_features.h>

#include "spatial_index.h"

#include <algorithm>
#include <cstddef>
#include <mutex>
#include <optional>
#include <span>
#include <utility>
#include <vector>

namespace chargefw::features {
namespace {

// The pair that std::map insertion in atom order would find first: the first atom that has an
// earlier atom at its position, and the first atom at that position. Returning there keeps many
// atoms at one position from costing a search each. Requires finite coordinates.
auto find_coincident_atoms(const ConformerFeatures& geometry,
                           const detail::ConformerSpatialIndex& index)
    -> std::optional<std::pair<std::size_t, std::size_t>> {
    auto neighbors = std::vector<std::size_t>{};
    for (std::size_t atom = 0; atom < geometry.molecule().atom_count(); ++atom) {
        const auto& position = geometry.position(atom);
        index.neighbor_indices_within(position, 0.0, neighbors);
        auto first = atom;
        for (const auto other : neighbors) {
            if (other < first && geometry.position(other) == position) {
                first = other;
            }
        }
        if (first != atom) {
            return std::pair{first, atom};
        }
    }
    return std::nullopt;
}

// Finds identical positions and atoms that are too close with one radius search per atom. Each
// atom searches within close_atom_factor times twice its own covalent radius, so a pair is always
// found from its atom with the larger radius and counted only from there. An element without a
// covalent radius counts as radius 0. At the first identical positions, the search stops and
// find_coincident_atoms() names the pair. Requires finite coordinates.
auto add_spatial_defects(const ConformerFeatures& geometry, ConformerGeometryDefects& defects)
    -> void {
    const auto& table = core::periodic_table();
    auto radii = std::vector<double>{};
    radii.reserve(geometry.molecule().atom_count());
    for (const auto& atom : geometry.molecule().atoms()) {
        radii.push_back(table.element(atom.atomic_number()).covalent_radius);
    }

    const auto index = detail::ConformerSpatialIndex{geometry};
    auto close = CloseAtoms{};
    auto closest_ratio = 0.0;
    auto neighbors = std::vector<std::size_t>{};
    for (std::size_t atom = 0; atom < radii.size(); ++atom) {
        index.neighbor_indices_within(geometry.position(atom),
                                      2.0 * close_atom_factor * radii[atom], neighbors);
        for (const auto other : neighbors) {
            if (other == atom) {
                continue;
            }
            if (geometry.position(other) == geometry.position(atom)) {
                defects.coincident_atoms = find_coincident_atoms(geometry, index);
                return;
            }
            if (radii[other] > radii[atom] || (radii[other] == radii[atom] && other < atom)) {
                continue;
            }
            const auto radius_sum = radii[atom] + radii[other];
            const auto distance = geometry.distance(atom, other);
            if (distance >= close_atom_factor * radius_sum) {
                continue;
            }
            const auto pair = std::pair{std::min(atom, other), std::max(atom, other)};
            const auto ratio = distance / radius_sum;
            if (close.pair_count == 0 || ratio < closest_ratio ||
                (ratio == closest_ratio && pair < close.closest)) {
                close.closest = pair;
                close.distance = distance;
                closest_ratio = ratio;
            }
            ++close.pair_count;
        }
    }
    if (close.pair_count > 0) {
        defects.close_atoms = close;
    }
}

} // namespace

auto PreparedMolecule::molecule() const noexcept -> const core::Molecule& {
    return topology_.molecule();
}

auto PreparedMolecule::topology() const noexcept -> const TopologyFeatures& {
    return topology_;
}

auto PreparedMolecule::geometry_defects() const -> std::span<const ConformerGeometryDefects> {
    // Built aside and stored only when complete: if this throws, the next call starts again.
    std::call_once(geometry_->found, [this] {
        const auto& source = molecule();
        auto conformers = std::vector<ConformerGeometryDefects>{};
        conformers.reserve(source.conformer_count());
        for (std::size_t conformer = 0; conformer < source.conformer_count(); ++conformer) {
            const auto geometry = ConformerFeatures{source, conformer};
            auto defects =
                ConformerGeometryDefects{.nonfinite_atom = geometry.first_nonfinite_atom_index()};
            if (!defects.nonfinite_atom.has_value()) {
                add_spatial_defects(geometry, defects);
            }
            conformers.push_back(defects);
        }
        geometry_->conformers = std::move(conformers);
    });
    return geometry_->conformers;
}

} // namespace chargefw::features
