#pragma once

#include "calculation/fixed_charge_partition.h"

#include <chargefw/calculation/assessment.h>
#include <chargefw/calculation/input_issues.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method_applicability.h>
#include <chargefw/parameters/models/parameter_set.h>

#include <memory>
#include <utility>
#include <vector>

namespace chargefw::calculation::detail {

// What an assessment found, shared by the assessment and all its plans so that a plan stays valid
// as long as anyone holds it. Never moved: the prepared molecules refer to the molecules here, and
// the plans to the applicable methods. Filled during assessment and only read afterwards.
struct AssessmentState {
    AssessmentState(core::MoleculeCollection supplied_molecules,
                    std::vector<parameters::ParameterSet> supplied_parameter_sets,
                    std::unique_ptr<FixedChargePartition> partition)
        : molecules{std::move(supplied_molecules)},
          parameter_sets{std::move(supplied_parameter_sets)},
          fixed_charge_partition{std::move(partition)},
          prepared_molecules{fixed_charge_partition ? fixed_charge_partition->active_molecules
                                                    : molecules} {}

    AssessmentState(const AssessmentState&) = delete;
    auto operator=(const AssessmentState&) -> AssessmentState& = delete;
    AssessmentState(AssessmentState&&) = delete;
    auto operator=(AssessmentState&&) -> AssessmentState& = delete;
    ~AssessmentState() = default;

    core::MoleculeCollection molecules;
    std::vector<parameters::ParameterSet> parameter_sets;
    std::unique_ptr<FixedChargePartition> fixed_charge_partition;
    // With fixed ions, the active parts of the molecules.
    features::PreparedMoleculeCollection prepared_molecules;
    methods::ApplicabilityResult applicability;
    std::vector<Rejection> rejections;
    // Every input issue of the molecules, whichever plan is executed.
    std::vector<InputIssue> input_warnings;
    double applicability_seconds = 0.0;
};

} // namespace chargefw::calculation::detail
