#pragma once

#include <chargefw/adapters/charge_result.h>

#include <iosfwd>

namespace chargefw::adapters::native::mol2_output {

// Writes fresh MOL2 records generated from the exact inputs retained by a calculation result.
class Mol2Writer {
  public:
    explicit Mol2Writer(std::ostream& output);

    auto write(const ChargeCalculationResult& result) const -> void;

  private:
    std::ostream* output_;
};

} // namespace chargefw::adapters::native::mol2_output
