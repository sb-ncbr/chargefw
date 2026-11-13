#pragma once

#include "calculation/fixed_charge_partition.h"
#include "methods/method_applicability.h"
#include "parameters/classification/classification_checks.h"

#include <chargefw/core/molecule.h>
#include <chargefw/features/conformer_features.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/methods/calculation_input.h>
#include <chargefw/methods/method.h>
#include <chargefw/parameters/classification/parameter_classification.h>
#include <chargefw/parameters/models/parameter_view.h>

#include <cstddef>
#include <memory>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace chargefw::calculation::detail {

// Builds the input of one full-execution target, a molecule or one of its conformers, and calls
// `operation(input)` while the objects the input borrows are alive: the conformer geometry, the
// fixed sources and target charge from `partition_target` (or the molecule's formal charge), and
// the parameter view of the checked `classification`. Shared by calculation and derivatives.
template <typename Operation>
auto with_target_input(const methods::ApplicableMethod& selected,
                       const features::PreparedMolecule& molecule,
                       const parameters::ParameterClassification* classification,
                       const std::optional<std::size_t> conformer_index,
                       const FixedChargePartitionTarget* partition_target, Operation&& operation)
    -> decltype(auto) {
    auto geometry = std::optional<features::ConformerFeatures>{};
    if (selected.method->requirements().coordinates) {
        if (!conformer_index.has_value()) {
            throw std::invalid_argument{"geometry-dependent method requires a conformer index"};
        }
        geometry.emplace(molecule.molecule(), *conformer_index);
    }

    auto fixed_sources = std::vector<methods::FixedPointSource>{};
    auto target_charge = core::total_formal_charge(molecule.molecule());
    if (partition_target != nullptr) {
        target_charge = partition_target->active_charge;
        fixed_sources = materialize_fixed_point_sources(*partition_target, conformer_index);
    }

    auto parameter_view = std::optional<parameters::ParameterView>{};
    if (classification != nullptr) {
        parameters::validate_parameter_classification(molecule.molecule(), *selected.parameter_set,
                                                      *classification);
        parameter_view.emplace(*selected.parameter_set, *classification);
    }

    const methods::CalculationInput input{molecule,
                                          selected.method_options,
                                          target_charge,
                                          geometry ? std::addressof(*geometry) : nullptr,
                                          parameter_view ? std::addressof(*parameter_view)
                                                         : nullptr,
                                          fixed_sources};
    return std::forward<Operation>(operation)(input);
}

} // namespace chargefw::calculation::detail
