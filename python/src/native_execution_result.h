#pragma once

#include <chargefw/io/charge_result.h>

#include <utility>

namespace chargefw::python {

class NativeExecutionResult {
  public:
    explicit NativeExecutionResult(io::ChargeCalculationResult result)
        : result_{std::move(result)} {}

    [[nodiscard]] auto result() const noexcept -> const io::ChargeCalculationResult& {
        return result_;
    }

  private:
    io::ChargeCalculationResult result_;
};

} // namespace chargefw::python
