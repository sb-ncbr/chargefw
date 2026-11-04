#pragma once

#include <chargefw/core/molecule.h>
#include <chargefw/core/position.h>

#include <cstddef>
#include <optional>

namespace chargefw::features {

class ConformerFeatures {
  public:
    explicit ConformerFeatures(const core::Molecule& molecule, std::size_t conformer_index = 0);

    explicit ConformerFeatures(core::Molecule&&) = delete;
    explicit ConformerFeatures(const core::Molecule&&) = delete;
    explicit ConformerFeatures(core::Molecule&&, std::size_t conformer_index) = delete;
    explicit ConformerFeatures(const core::Molecule&&, std::size_t conformer_index) = delete;

    [[nodiscard]] auto molecule() const noexcept -> const core::Molecule&;
    [[nodiscard]] auto conformer_index() const noexcept -> std::size_t;

    [[nodiscard]] auto position(std::size_t atom_index) const -> const core::Position&;

    [[nodiscard]] auto squared_distance(std::size_t first_atom_index,
                                        std::size_t second_atom_index) const -> double;
    [[nodiscard]] auto distance(std::size_t first_atom_index, std::size_t second_atom_index) const
        -> double;

    [[nodiscard]] auto first_nonfinite_atom_index() const noexcept -> std::optional<std::size_t>;

  private:
    auto validate_conformer_index() const -> void;
    auto validate_atom_index(std::size_t atom_index) const -> void;

    const core::Molecule* molecule_;
    std::size_t conformer_index_;
};

} // namespace chargefw::features
