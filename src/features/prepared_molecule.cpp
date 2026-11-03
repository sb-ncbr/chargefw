#include <chargefw/core/molecule.h>
#include <chargefw/features/conformer_features.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/features/topology_features.h>

#include <cstddef>
#include <mutex>
#include <span>
#include <utility>
#include <vector>

namespace chargefw::features {

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
                defects.coincident_atoms = geometry.coincident_atom_indices();
            }
            conformers.push_back(defects);
        }
        geometry_->conformers = std::move(conformers);
    });
    return geometry_->conformers;
}

} // namespace chargefw::features