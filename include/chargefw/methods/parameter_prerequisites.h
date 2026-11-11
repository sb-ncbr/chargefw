#pragma once

#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method_prerequisites.h>
#include <chargefw/parameters/classification/classification_result.h>
#include <chargefw/parameters/classification/parameter_classification.h>
#include <chargefw/parameters/models/parameter_set.h>

#include <vector>

namespace chargefw::methods {

class Method;

struct CollectionParameterPrerequisiteResult {
    std::vector<PrerequisiteIssue> issues;
    std::vector<parameters::ParameterClassification> classifications;

    [[nodiscard]] explicit operator bool() const noexcept {
        return issues.empty();
    }
};

[[nodiscard]] auto
check_parameter_prerequisites(const Method& method,
                              const features::PreparedMoleculeCollection& molecules,
                              const parameters::ParameterSet& parameter_set,
                              const parameters::ClassificationOptions& classification_options = {})
    -> CollectionParameterPrerequisiteResult;

} // namespace chargefw::methods
