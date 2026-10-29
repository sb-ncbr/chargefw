#pragma once

#include <chargefw/adapters/charge_result.h>

#include <iosfwd>
#include <span>
#include <string_view>

namespace chargefw::adapters::gemmi::mmcif_output {

// Writes a fresh mmCIF representation of calculation results.
class MmcifWriter {
  public:
    explicit MmcifWriter(std::ostream& output);

    auto write(const ChargeCalculationResult& result, std::string_view generator_version = {}) const
        -> void;
    // Writes several charge sets for the same molecules; charge set k gets metadata id k + 1.
    // Results must have equal inputs and distinct method and parameter set pairs.
    auto write(std::span<const ChargeCalculationResult* const> results,
               std::string_view generator_version = {}) const -> void;

  private:
    std::ostream* output_;
};

} // namespace chargefw::adapters::gemmi::mmcif_output
