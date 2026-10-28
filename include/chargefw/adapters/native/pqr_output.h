#pragma once

#include <chargefw/adapters/charge_result.h>

#include <iosfwd>
#include <string_view>

namespace chargefw::adapters::native::pqr_output {

// Writes PQR for one structure imported from PDB or mmCIF: source atom labels, coordinates,
// charges, and the van der Waals radius of each atom's element.
class PqrWriter {
  public:
    explicit PqrWriter(std::ostream& output);

    auto write(const ChargeCalculationResult& result, std::string_view generator_version = {}) const
        -> void;

  private:
    std::ostream* output_;
};

} // namespace chargefw::adapters::native::pqr_output
