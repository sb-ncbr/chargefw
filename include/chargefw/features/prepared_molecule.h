#pragma once

#include <chargefw/core/molecule.h>
#include <chargefw/features/topology_features.h>

#include <cstddef>
#include <memory>
#include <mutex>
#include <optional>
#include <span>
#include <utility>
#include <vector>

namespace chargefw::features {

// Atoms are too close when their distance is below this fraction of the sum of their covalent
// radii. Bonds in clean structures stay above 0.78 of it.
inline constexpr double close_atom_factor = 0.7;

// The pairs of atoms of one conformer that are too close, but not at identical positions.
struct CloseAtoms {
    // The pair with the smallest distance relative to the sum of its covalent radii.
    std::pair<std::size_t, std::size_t> closest;
    double distance = 0.0; // Å
    std::size_t pair_count = 0;
};

// Defects of one conformer's geometry that matter to methods using coordinates.
struct ConformerGeometryDefects {
    // The first atom with a non-finite coordinate. Methods using coordinates cannot calculate.
    std::optional<std::size_t> nonfinite_atom;
    // Two atoms at identical positions, the later atom being the first one in atom order that
    // has an earlier atom at its position. Methods using coordinates cannot calculate. Looked for
    // only when every coordinate is finite.
    std::optional<std::pair<std::size_t, std::size_t>> coincident_atoms = std::nullopt;
    // Methods using coordinates can calculate, but their charges may be wrong. Looked for only
    // when there is no defect above.
    std::optional<CloseAtoms> close_atoms = std::nullopt;
};

class PreparedMolecule {
  public:
    explicit PreparedMolecule(const core::Molecule& molecule)
        : topology_{molecule}, geometry_{std::make_shared<GeometryCache>()} {}

    PreparedMolecule(core::Molecule&&) = delete;
    PreparedMolecule(const core::Molecule&&) = delete;

    // Copies refer to the same molecule, so they share its geometry defects. Moving shares them too
    // instead of taking them, so that geometry queries remain valid on the moved-from object.
    PreparedMolecule(const PreparedMolecule&) = default;
    auto operator=(const PreparedMolecule&) -> PreparedMolecule& = default;
    PreparedMolecule(PreparedMolecule&& other) noexcept
        // NOLINTNEXTLINE(performance-move-constructor-init): the cache is shared, not taken.
        : topology_{std::move(other.topology_)}, geometry_{other.geometry_} {}
    auto operator=(PreparedMolecule&& other) noexcept -> PreparedMolecule& {
        if (this != &other) {
            topology_ = std::move(other.topology_);
            geometry_ = other.geometry_;
        }
        return *this;
    }
    ~PreparedMolecule() = default;

    [[nodiscard]] auto molecule() const noexcept -> const core::Molecule&;

    [[nodiscard]] auto topology() const noexcept -> const TopologyFeatures&;

    // The geometry defects of each conformer, found on first use and then reused; safe to call
    // from several threads. Molecules that no method using coordinates looks at never find them.
    [[nodiscard]] auto geometry_defects() const -> std::span<const ConformerGeometryDefects>;

  private:
    struct GeometryCache {
        std::once_flag found;
        std::vector<ConformerGeometryDefects> conformers;
    };

    TopologyFeatures topology_;
    // Behind a pointer because std::once_flag cannot move, and prepared molecules live in a vector.
    // Never null.
    std::shared_ptr<GeometryCache> geometry_;
};

} // namespace chargefw::features
