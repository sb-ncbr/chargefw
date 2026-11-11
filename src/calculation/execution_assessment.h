#pragma once

#include <chargefw/calculation/execution_issues.h>
#include <chargefw/calculation/execution_policy.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method.h>

#include <cstdint>
#include <vector>

namespace chargefw::calculation::detail {

enum class ExecutionAvailability : std::uint8_t {
    available,
    available_with_warning,
    unsupported,
};

struct ExecutionAssessment {
    ExecutionAvailability availability;
    std::vector<ExecutionIssue> issues;
};

// Whether `method` can calculate `molecules` in `mode`, and what the resource thresholds say about
// it. Depends only on what the method declares, not on its parameters.
[[nodiscard]] auto assess_execution(const methods::Method& method,
                                    const features::PreparedMoleculeCollection& molecules,
                                    const ResourcePolicy& resource_policy, ExecutionMode mode)
    -> ExecutionAssessment;

} // namespace chargefw::calculation::detail
