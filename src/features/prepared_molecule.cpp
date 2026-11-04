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

// Counts the atoms that are too close with one radius search per atom. Each atom searches within
// close_atom_factor times twice its own covalent radius, so a pair is always found from its atom
// with the larger radius and counted only from there. An element without a covalent radius counts
// as radius 0. Atoms at identical positions count as a pair at distance 0, unless
// stop_at_identical_positions, which returns std::nullopt at the first of them instead. Requires
// finite coordinates.
auto count_close_atoms(const ConformerFeatures& geometry,
                       const detail::ConformerSpatialIndex& index,
                       const bool stop_at_identical_positions) -> std::optional<CloseAtoms> {
    const auto& table = core::periodic_table();
    auto radii = std::vector<double>{};
    radii.reserve(geometry.molecule().atom_count());
    for (const auto& atom : geometry.molecule().atoms()) {
        radii.push_back(table.element(atom.atomic_number()).covalent_radius);
    }

    auto close = CloseAtoms{};
    auto closest_ratio = 0.0;
    auto neighbors = std::vector<std::size_t>{};
    // Allocated at the first pair; marks rather than a list keep many pairs of one atom from
    // growing it beyond the atom count.
    auto in_pair = std::vector<bool>{};
    for (std::size_t atom = 0; atom < radii.size(); ++atom) {
        index.neighbor_indices_within(geometry.position(atom),
                                      2.0 * close_atom_factor * radii[atom], neighbors);
        for (const auto other : neighbors) {
            if (other == atom) {
                continue;
            }
            if (stop_at_identical_positions &&
                geometry.position(other) == geometry.position(atom)) {
                return std::nullopt;
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
            in_pair.resize(radii.size());
            in_pair[atom] = true;
            in_pair[other] = true;
        }
    }
    for (std::size_t atom = 0; atom < in_pair.size(); ++atom) {
        if (in_pair[atom]) {
            close.atoms.push_back(atom);
        }
    }
    return close;
}

auto find_geometry_defects(const ConformerFeatures& geometry) -> ConformerGeometryDefects {
    auto defects =
        ConformerGeometryDefects{.nonfinite_atom = geometry.first_nonfinite_atom_index()};
    if (defects.nonfinite_atom.has_value()) {
        return defects;
    }
    const auto index = detail::ConformerSpatialIndex{geometry};
    const auto close = count_close_atoms(geometry, index, true);
    if (!close.has_value()) {
        defects.coincident_atoms = find_coincident_atoms(geometry, index);
    } else if (close->pair_count > 0) {
        defects.close_atoms = close;
    }
    return defects;
}

} // namespace

auto find_close_atoms(const ConformerFeatures& geometry) -> std::optional<CloseAtoms> {
    const auto index = detail::ConformerSpatialIndex{geometry};
    const auto close = count_close_atoms(geometry, index, false);
    return close->pair_count > 0 ? close : std::nullopt;
}

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
            conformers.push_back(find_geometry_defects(ConformerFeatures{source, conformer}));
        }
        geometry_->conformers = std::move(conformers);
    });
    return geometry_->conformers;
}

} // namespace chargefw::features
