#pragma once

#include "calculation/assessment_state.h"
#include "calculation/fixed_charge_partition.h"
#include "methods/method_applicability.h"

#include <chargefw/calculation/assessment.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/parameters/classification/parameter_classification.h>
#include <chargefw/parameters/models/parameter_set.h>

#include <span>

namespace chargefw::calculation::detail {

// What other parts of the library may read from a plan's assessment, so that its public header
// names none of them. Read-only.
struct AssessmentAccess {
    // The molecules of the assessment, as supplied.
    [[nodiscard]] static auto source_molecules(const ExecutionPlan& plan)
        -> const core::MoleculeCollection& {
        return plan.state_->molecules;
    }

    // The prepared molecules the plan calculates: with fixed ions, the active parts.
    [[nodiscard]] static auto molecules(const ExecutionPlan& plan)
        -> const features::PreparedMoleculeCollection& {
        return plan.state_->prepared_molecules;
    }

    // The fixed-ion partition, or nullptr without fixed ions.
    [[nodiscard]] static auto partition(const ExecutionPlan& plan) -> const FixedChargePartition* {
        return plan.state_->fixed_charge_partition.get();
    }

    // The parameter entries matched to the atoms and bonds the plan calculates, one classification
    // per calculated molecule (ExecutionPlan::classifications()).
    [[nodiscard]] static auto classifications(const ExecutionPlan& plan)
        -> std::span<const parameters::ParameterClassification> {
        return plan.classifications();
    }

    // The method of `plan` with the values of `candidate`, checked as calculate(plan, candidate)
    // checks them.
    [[nodiscard]] static auto with_values(const ExecutionPlan& plan,
                                          const parameters::ParameterSet& candidate)
        -> methods::ApplicableMethod {
        return *plan.with_values(&candidate);
    }
};

} // namespace chargefw::calculation::detail
