#pragma once

#include <chargefw/adapters/charge_result.h>

#include <iosfwd>
#include <optional>

namespace chargefw::adapters::native::json_output {

// Writes a version 1.0 ChargeFW calculation-result JSON document.
class JsonWriter {
  public:
    explicit JsonWriter(std::ostream& output);

    auto write(const ChargeCalculationResult& result,
               const std::optional<ExecutionMetrics>& execution_metrics = std::nullopt) const
        -> void;

  private:
    std::ostream* output_;
};

} // namespace chargefw::adapters::native::json_output
