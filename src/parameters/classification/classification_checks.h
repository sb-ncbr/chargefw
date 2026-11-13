#pragma once

#include <chargefw/core/molecule.h>
#include <chargefw/parameters/classification/parameter_classification.h>
#include <chargefw/parameters/models/parameter_set.h>

namespace chargefw::parameters {

// Throws std::invalid_argument unless the classification has one entry index per atom when the set
// has atom entries and none otherwise, the same for bonds, and every index names an entry of the
// set.
auto validate_parameter_classification(const core::Molecule& molecule,
                                       const ParameterSet& parameters,
                                       const ParameterClassification& classification) -> void;

} // namespace chargefw::parameters
