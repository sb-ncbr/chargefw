#pragma once

#include <chargefw/adapters/charge_result.h>

#include <iosfwd>
#include <span>

#include <gemmi/cifdoc.hpp>

namespace chargefw::adapters::gemmi::mmcif_output {

auto write_attached(std::ostream& output, std::span<const ChargeCalculationResult* const> results,
                    const ::gemmi::cif::Document& source, bool overwrite) -> void;

} // namespace chargefw::adapters::gemmi::mmcif_output
