#pragma once

#include <chargefw/io/charge_result.h>

#include <cstdint>
#include <iosfwd>
#include <span>
#include <string_view>

namespace chargefw::io {

// What annotation does with partial charge categories that the source already has.
enum class ExistingCharges : std::uint8_t {
    error,
    replace,
    // Keeps the charge sets of the source and gives the new ones the next IDs. The source's sets
    // must follow dictionary 1.1 and differ from the new ones in method or parameter set.
    append,
};

[[nodiscard]] auto existing_charges_from_string(std::string_view value) -> ExistingCharges;

// Writes calculation results as mmCIF: write() generates the structure from the calculated
// molecules; annotate() adds the charges to the structure the molecules were read from.
class MmcifWriter {
  public:
    explicit MmcifWriter(std::ostream& output);

    // Writes the mmCIF input the molecules were read from with their charges added. Only the
    // calculated atom sites get charge rows, which refer to them by _atom_site.id; every other
    // category, block, and atom site of the input is kept. The results must be for the same
    // molecules, read from one mmCIF input.
    auto annotate(std::span<const ChargeCalculationResult* const> results,
                  ExistingCharges existing = ExistingCharges::error) const -> void;

    auto write(const ChargeCalculationResult& result) const -> void;
    // Writes several charge sets for the same molecules; charge set k gets metadata id k + 1.
    // Results must have equal inputs and distinct method and parameter set pairs.
    auto write(std::span<const ChargeCalculationResult* const> results) const -> void;

  private:
    std::ostream* output_;
};

} // namespace chargefw::io
