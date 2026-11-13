#pragma once

#include "methods/method_applicability.h"

#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/parameters/classification/classification_result.h>
#include <chargefw/parameters/models/parameter_set.h>

#include <span>
#include <string>

namespace chargefw::methods {

class Method;

struct ApplicabilityRequest {
    // NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
    const features::PreparedMoleculeCollection& molecules;
    std::span<const Method* const> methods;
    std::span<const parameters::ParameterSet> parameter_sets;
    // NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)
    parameters::ClassificationOptions classification_options{};
    // Applied over each method's defaults. Options come only with a named method, which is then the
    // only one in `methods`; with automatic selection they are empty.
    MethodOptions method_options;
};

[[nodiscard]] auto find_applicable_methods(const ApplicabilityRequest& request)
    -> ApplicabilityResult;

} // namespace chargefw::methods
