#include <chargefw/core/position.h>
#include <chargefw/features/conformer_features.h>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace chargefw::features {

ConformerFeatures::ConformerFeatures(const core::Molecule& molecule,
                                     const std::size_t conformer_index)
    : molecule_{&molecule}, conformer_index_{conformer_index} {
    validate_conformer_index();
}

auto ConformerFeatures::molecule() const noexcept -> const core::Molecule& {
    return *molecule_;
}

auto ConformerFeatures::conformer_index() const noexcept -> std::size_t {
    return conformer_index_;
}

auto ConformerFeatures::position(const std::size_t atom_index) const -> const core::Position& {
    validate_atom_index(atom_index);
    return molecule_->conformer(conformer_index_)[atom_index];
}

auto ConformerFeatures::squared_distance(const std::size_t first_atom_index,
                                         const std::size_t second_atom_index) const -> double {
    return core::squared_distance(position(first_atom_index), position(second_atom_index));
}

auto ConformerFeatures::distance(const std::size_t first_atom_index,
                                 const std::size_t second_atom_index) const -> double {
    return core::distance(position(first_atom_index), position(second_atom_index));
}

auto ConformerFeatures::first_nonfinite_atom_index() const noexcept -> std::optional<std::size_t> {
    const auto positions = molecule_->conformer(conformer_index_).positions();
    const auto found = std::ranges::find_if(positions, [](const core::Position& position) -> bool {
        return !std::isfinite(position.x) || !std::isfinite(position.y) ||
               !std::isfinite(position.z);
    });

    if (found == positions.end()) {
        return std::nullopt;
    }
    return static_cast<std::size_t>(found - positions.begin());
}

auto ConformerFeatures::validate_conformer_index() const -> void {
    if (conformer_index_ >= molecule_->conformer_count()) {
        throw std::out_of_range{"conformer index is outside the molecule"};
    }
}

auto ConformerFeatures::validate_atom_index(const std::size_t atom_index) const -> void {
    if (atom_index >= molecule_->atom_count()) {
        throw std::out_of_range{"atom index is outside the molecule"};
    }
}

} // namespace chargefw::features
