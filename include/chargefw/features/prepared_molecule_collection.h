#pragma once

#include <chargefw/core/molecule_collection.h>
#include <chargefw/features/prepared_molecule.h>

#include <cstddef>
#include <span>
#include <vector>

namespace chargefw::features {

class PreparedMoleculeCollection {
  public:
    explicit PreparedMoleculeCollection(const core::MoleculeCollection& collection);

    explicit PreparedMoleculeCollection(core::MoleculeCollection&&) = delete;
    explicit PreparedMoleculeCollection(const core::MoleculeCollection&&) = delete;

    // Moving keeps each prepared molecule in place; copying would need to copy them.
    PreparedMoleculeCollection(const PreparedMoleculeCollection&) = delete;
    auto operator=(const PreparedMoleculeCollection&) -> PreparedMoleculeCollection& = delete;
    PreparedMoleculeCollection(PreparedMoleculeCollection&&) noexcept = default;
    auto operator=(PreparedMoleculeCollection&&) noexcept -> PreparedMoleculeCollection& = default;
    ~PreparedMoleculeCollection() = default;

    [[nodiscard]] auto molecules() const noexcept -> std::span<const PreparedMolecule>;

    [[nodiscard]] auto size() const noexcept -> std::size_t;
    [[nodiscard]] auto empty() const noexcept -> bool;

    [[nodiscard]] auto operator[](std::size_t index) const noexcept -> const PreparedMolecule&;
    [[nodiscard]] auto at(std::size_t index) const -> const PreparedMolecule&;

  private:
    std::vector<PreparedMolecule> molecules_;
};
} // namespace chargefw::features