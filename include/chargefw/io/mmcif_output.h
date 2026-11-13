#pragma once

#include <chargefw/io/charge_result.h>

#include <iosfwd>
#include <span>

namespace chargefw::io {

// Writes a fresh mmCIF representation of calculation results.
class MmcifWriter {
  public:
    explicit MmcifWriter(std::ostream& output);

    auto write(const ChargeCalculationResult& result) const -> void;
    // Writes several charge sets for the same molecules; charge set k gets metadata id k + 1.
    // Results must have equal inputs and distinct method and parameter set pairs.
    auto write(std::span<const ChargeCalculationResult* const> results) const -> void;

  private:
    std::ostream* output_;
};

} // namespace chargefw::io
