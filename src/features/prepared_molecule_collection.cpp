#include <chargefw/features/prepared_molecule.h>
#include <chargefw/features/prepared_molecule_collection.h>

namespace chargefw::features {

// Built from the range in one step, which constructs each element in place; adding them one by one
// would need prepared molecules to be movable.
PreparedMoleculeCollection::PreparedMoleculeCollection(const core::MoleculeCollection& collection)
    : molecules_(collection.molecules().begin(), collection.molecules().end()) {}

auto PreparedMoleculeCollection::molecules() const noexcept -> std::span<const PreparedMolecule> {
    return molecules_;
}

auto PreparedMoleculeCollection::size() const noexcept -> std::size_t {
    return molecules_.size();
}

auto PreparedMoleculeCollection::empty() const noexcept -> bool {
    return molecules_.empty();
}

auto PreparedMoleculeCollection::operator[](std::size_t index) const noexcept
    -> const PreparedMolecule& {
    return molecules_[index];
}

auto PreparedMoleculeCollection::at(std::size_t index) const -> const PreparedMolecule& {
    return molecules_.at(index);
}

} // namespace chargefw::features