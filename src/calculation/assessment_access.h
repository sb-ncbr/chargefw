#pragma once

#include "calculation/fixed_charge_partition.h"

#include <chargefw/calculation/assessment.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method_applicability.h>
#include <chargefw/parameters/classification/parameter_classification.h>
#include <chargefw/parameters/models/parameter_set.h>

#include <span>

namespace chargefw::calculation::detail {

// What other parts of the library may read from an assessment, so that its public header names
// none of them. Read-only; a function that takes a plan checks that it belongs to the assessment.
struct AssessmentAccess {
    // The prepared molecules the plans calculate: with fixed ions, the active parts.
    [[nodiscard]] static auto molecules(const AssessmentResult& assessment)
        -> const features::PreparedMoleculeCollection& {
        return assessment.prepared_molecules();
    }

    // The fixed-ion partition, or nullptr without fixed ions.
    [[nodiscard]] static auto partition(const AssessmentResult& assessment)
        -> const FixedChargePartition* {
        return assessment.fixed_charge_partition_.get();
    }

    // The parameter entries matched to the atoms and bonds `plan` calculates, one classification
    // per calculated molecule (ExecutionPlan::classifications()).
    [[nodiscard]] static auto classifications(const AssessmentResult& assessment,
                                              const ExecutionPlan& plan)
        -> std::span<const parameters::ParameterClassification> {
        assessment.check_plan(plan);
        return plan.classifications();
    }

    // The method of `plan` with the values of `candidate`, checked as calculate(assessment, plan,
    // candidate) checks them.
    [[nodiscard]] static auto with_values(const AssessmentResult& assessment,
                                          const ExecutionPlan& plan,
                                          const parameters::ParameterSet& candidate)
        -> methods::ApplicableMethod {
        assessment.check_plan(plan);
        return *AssessmentResult::replace_values(plan, &candidate);
    }
};

} // namespace chargefw::calculation::detail
