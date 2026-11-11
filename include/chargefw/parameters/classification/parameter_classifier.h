#pragma once

#include <chargefw/features/topology_features.h>
#include <chargefw/parameters/classification/classification_result.h>
#include <chargefw/parameters/classification/parameter_classification.h>
#include <chargefw/parameters/models/parameter_set.h>

namespace chargefw::parameters {

[[nodiscard]] auto try_classify_parameters(const features::TopologyFeatures& topology,
                                           const ParameterSet& parameters,
                                           const ClassificationOptions& options = {})
    -> ClassificationResult;

[[nodiscard]] auto classify_parameters(const features::TopologyFeatures& topology,
                                       const ParameterSet& parameters,
                                       const ClassificationOptions& options = {})
    -> ParameterClassification;

} // namespace chargefw::parameters
