#pragma once

#include "features/topology_helpers.h"

#include <chargefw/core/molecule.h>
#include <chargefw/features/topology_features.h>

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace chargefw::methods::detail {

// Shifts reference charges, one per atom, uniformly within each connected component so that the
// component sums to the formal charge of its atoms. Charge cannot move between components in
// methods that transfer it only along bonds, so each component keeps this sum. A target charge that
// differs from the total formal charge is then spread over all atoms; for a single component this
// is the uniform shift to the target charge.
[[nodiscard]] inline auto component_normalized_charges(const core::Molecule& molecule,
                                                       const features::TopologyFeatures& topology,
                                                       std::vector<double> charges,
                                                       const double target_charge)
    -> std::vector<double> {
    if (charges.size() != molecule.atom_count()) {
        throw std::logic_error{"reference charge count must match atom count"};
    }
    if (charges.empty()) {
        return charges;
    }
    auto formal_total = 0.0;
    for (const auto& component : features::connected_components(topology.adjacency())) {
        auto reference = 0.0;
        auto formal = 0.0;
        for (const auto atom_index : component) {
            reference += charges[atom_index];
            formal += molecule.atom(atom_index).formal_charge();
        }
        const auto shift = (formal - reference) / static_cast<double>(component.size());
        for (const auto atom_index : component) {
            charges[atom_index] += shift;
        }
        formal_total += formal;
    }
    const auto remainder = (target_charge - formal_total) / static_cast<double>(charges.size());
    if (remainder != 0.0) {
        for (auto& charge : charges) {
            charge += remainder;
        }
    }
    return charges;
}

// The linear part of component_normalized_charges(): subtracts each connected component's mean.
// The normalization is this projection plus a constant, so a change of the reference charges
// changes the normalized charges by the projected change. The projection is symmetric, so it also
// maps a gradient with respect to the normalized charges to one with respect to the references.
// `components` are the atom indices of each connected component, as
// features::connected_components() returns them; a caller that outlives the topology keeps this
// list.
inline auto remove_component_means(const std::vector<std::vector<std::size_t>>& components,
                                   std::vector<double>& values) -> void {
    for (const auto& component : components) {
        auto sum = 0.0;
        for (const auto atom_index : component) {
            sum += values[atom_index];
        }
        const auto mean = sum / static_cast<double>(component.size());
        for (const auto atom_index : component) {
            values[atom_index] -= mean;
        }
    }
}

} // namespace chargefw::methods::detail
