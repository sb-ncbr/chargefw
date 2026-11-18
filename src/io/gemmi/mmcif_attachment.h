#pragma once

#include <chargefw/io/charge_result.h>
#include <chargefw/io/mmcif_output.h>

#include <iosfwd>
#include <span>

#include <gemmi/cifdoc.hpp>

namespace chargefw::io::gemmi::mmcif_output {

// Adds the charges to the document and writes it. Every block with coordinates must belong to a
// calculated molecule.
auto write_attached(std::ostream& output, std::span<const ChargeCalculationResult* const> results,
                    ::gemmi::cif::Document document, ExistingCharges existing) -> void;

} // namespace chargefw::io::gemmi::mmcif_output
