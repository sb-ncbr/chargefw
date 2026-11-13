#pragma once

#include <chargefw/io/charge_result.h>

#include <iosfwd>

namespace chargefw::io {

// Writes PQR for one structure imported from PDB or mmCIF: source atom labels, coordinates,
// charges, and the van der Waals radius of each atom's element.
class PqrWriter {
  public:
    explicit PqrWriter(std::ostream& output);

    auto write(const ChargeCalculationResult& result) const -> void;

  private:
    std::ostream* output_;
};

} // namespace chargefw::io
