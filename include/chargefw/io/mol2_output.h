#pragma once

#include <chargefw/io/charge_result.h>

#include <iosfwd>

namespace chargefw::io {

// Writes fresh MOL2 records generated from the exact inputs retained by a calculation result.
class Mol2Writer {
  public:
    explicit Mol2Writer(std::ostream& output);

    auto write(const ChargeCalculationResult& result) const -> void;

  private:
    std::ostream* output_;
};

} // namespace chargefw::io
