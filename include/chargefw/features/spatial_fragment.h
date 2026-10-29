#pragma once

#include <chargefw/core/molecule.h>
#include <chargefw/features/prepared_molecule.h>

#include <cstddef>
#include <memory>
#include <span>
#include <vector>

namespace chargefw::features {

class ConformerFeatures;
class SpatialFragmentBuilder;

class SpatialFragment {
  public:
    [[nodiscard]] auto molecule() const noexcept -> const core::Molecule&;
    [[nodiscard]] auto source_conformer_index() const noexcept -> std::size_t;
    [[nodiscard]] auto center_local_atom_index() const noexcept -> std::size_t;

    // Local atom order follows the source molecule's atom order.
    [[nodiscard]] auto local_to_source_atom_indices() const noexcept
        -> std::span<const std::size_t>;
    [[nodiscard]] auto local_to_source_bond_indices() const noexcept
        -> std::span<const std::size_t>;

  private:
    SpatialFragment(core::Molecule molecule, std::vector<std::size_t> local_to_source_atom_indices,
                    std::vector<std::size_t> local_to_source_bond_indices,
                    std::size_t center_local_atom_index, std::size_t source_conformer_index);

    core::Molecule molecule_;
    std::vector<std::size_t> local_to_source_atom_indices_;
    std::vector<std::size_t> local_to_source_bond_indices_;
    std::size_t center_local_atom_index_;
    std::size_t source_conformer_index_;

    friend class SpatialFragmentBuilder;
};

// Owns one KD-tree for a source conformer and builds its spatial fragments. It borrows the source
// molecule and geometry, which must outlive the builder.
class SpatialFragmentBuilder {
  public:
    SpatialFragmentBuilder(const PreparedMolecule& source, const ConformerFeatures& geometry);
    SpatialFragmentBuilder(const PreparedMolecule&&, const ConformerFeatures&) = delete;
    SpatialFragmentBuilder(const PreparedMolecule&, const ConformerFeatures&&) = delete;
    SpatialFragmentBuilder(const PreparedMolecule&&, const ConformerFeatures&&) = delete;
    ~SpatialFragmentBuilder();

    SpatialFragmentBuilder(const SpatialFragmentBuilder&) = delete;
    auto operator=(const SpatialFragmentBuilder&) -> SpatialFragmentBuilder& = delete;
    SpatialFragmentBuilder(SpatialFragmentBuilder&&) = delete;
    auto operator=(SpatialFragmentBuilder&&) -> SpatialFragmentBuilder& = delete;

    // Sorted source indices of atoms within the inclusive radius of the center atom.
    [[nodiscard]] auto atom_indices_within(std::size_t center_atom_index, double radius) const
        -> std::vector<std::size_t>;
    [[nodiscard]] auto build(std::size_t center_atom_index, double radius) const -> SpatialFragment;

  private:
    class SpatialIndex;

    const PreparedMolecule* source_;
    const ConformerFeatures* geometry_;
    std::unique_ptr<const SpatialIndex> spatial_index_;
};

} // namespace chargefw::features
