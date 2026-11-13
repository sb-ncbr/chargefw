#pragma once

#include "native_execution_result.h"

#include <chargefw/calculation/assessment.h>
#include <chargefw/io/charge_result.h>

#include <nanobind/nanobind.h>

#include <cstddef>
#include <memory>
#include <optional>
#include <utility>
#include <vector>

namespace chargefw::python {

// An assessment with what the Python results need, shared by its plans and by the objects that
// fit parameters on them (parameterization.cpp), which keep it alive.
class NativeAssessmentState {
  public:
    NativeAssessmentState(calculation::AssessmentResult assessment,
                          std::vector<io::ImportedMoleculeRecord> inputs,
                          io::RequestedCalculationProvenance requested,
                          const std::size_t max_threads)
        : assessment_{std::move(assessment)}, inputs_{std::move(inputs)},
          requested_{std::move(requested)}, max_threads_{max_threads} {}

    calculation::AssessmentResult assessment_;
    std::vector<io::ImportedMoleculeRecord> inputs_;
    io::RequestedCalculationProvenance requested_;
    std::size_t max_threads_ = 0;
};

// One plan of an assessment (calculation.cpp).
class NativePlan {
  public:
    NativePlan(std::shared_ptr<NativeAssessmentState> state, const std::size_t index)
        : state_{std::move(state)}, index_{index} {}

    [[nodiscard]] auto report() const -> nanobind::dict;
    [[nodiscard]] auto calculate(std::optional<std::size_t> max_threads,
                                 nanobind::object observer) const -> NativeExecutionResult;
    [[nodiscard]] auto estimate(std::optional<std::size_t> max_threads) const
        -> std::pair<double, double>;

    [[nodiscard]] auto state() const noexcept -> const std::shared_ptr<NativeAssessmentState>& {
        return state_;
    }
    [[nodiscard]] auto plan() const -> const calculation::ExecutionPlan& {
        return state_->assessment_.plans()[index_];
    }

  private:
    std::shared_ptr<NativeAssessmentState> state_;
    std::size_t index_ = 0;
};

} // namespace chargefw::python
