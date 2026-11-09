#pragma once

#include <chargefw/core/molecule.h>
#include <chargefw/features/topology_features.h>
#include <chargefw/parameters/models/parameter_key.h>

#include <cstddef>

// The keys the classifier gives single atoms and bonds under strict types, for code that builds
// parameter sets from molecules (parameterization::make_parameter_template()). Private to the
// library; the classifier stays the one place that decides types.
namespace chargefw::parameters::detail {

// Orders atom keys by atomic number, then classifier and type.
[[nodiscard]] auto atom_key_less(const AtomParameterKey& first, const AtomParameterKey& second)
    -> bool;

[[nodiscard]] auto atom_type_key(const core::Molecule& molecule,
                                 const features::TopologyFeatures& topology, std::size_t atom_index,
                                 AtomParameterClassificationKind classification)
    -> AtomParameterKey;

// A bond matches its key in either endpoint order, so the key puts the smaller atom key
// (atom_key_less()) first.
[[nodiscard]] auto bond_type_key(const core::Molecule& molecule,
                                 const features::TopologyFeatures& topology, std::size_t bond_index,
                                 AtomParameterClassificationKind atom_classification,
                                 BondParameterClassificationKind bond_classification)
    -> BondParameterKey;

} // namespace chargefw::parameters::detail
