#pragma once

#include <chargefw/core/molecule.h>

#include <cstddef>
#include <span>
#include <vector>

namespace chargefw::features {

class TopologyFeatures {
  public:
    explicit TopologyFeatures(const core::Molecule& molecule);

    explicit TopologyFeatures(core::Molecule&&) = delete;
    explicit TopologyFeatures(const core::Molecule&&) = delete;

    // Neither copied nor moved, like the prepared molecule that holds it.
    TopologyFeatures(const TopologyFeatures&) = delete;
    auto operator=(const TopologyFeatures&) -> TopologyFeatures& = delete;
    TopologyFeatures(TopologyFeatures&&) = delete;
    auto operator=(TopologyFeatures&&) -> TopologyFeatures& = delete;
    ~TopologyFeatures() = default;

    [[nodiscard]] auto molecule() const noexcept -> const core::Molecule&;

    [[nodiscard]] auto degree(std::size_t atom_index) const -> std::size_t;

    [[nodiscard]] auto neighbor_indices(std::size_t atom_index) const
        -> std::span<const std::size_t>;

    [[nodiscard]] auto incident_bond_indices(std::size_t atom_index) const
        -> std::span<const std::size_t>;

  private:
    const core::Molecule* molecule_;

    std::vector<std::vector<std::size_t>> neighbor_indices_;
    std::vector<std::vector<std::size_t>> incident_bond_indices_;
};

} // namespace chargefw::features
