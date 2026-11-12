#pragma once

#include <stdexcept>

namespace chargefw::calculation {

// The public calculation API converts this internal exception to ExecutionStatus::cancelled.
class CalculationCancelled : public std::runtime_error {
  public:
    CalculationCancelled() : std::runtime_error{"calculation was cancelled"} {}
};

} // namespace chargefw::calculation
