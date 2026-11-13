#pragma once

#include <chargefw/calculation/assessment.h>
#include <chargefw/io/molecule_record.h>

#include <span>
#include <string>
#include <vector>

namespace chargefw::io::detail {

[[nodiscard]] auto fixed_ion_names(bool common_only) -> std::vector<std::string>;

[[nodiscard]] auto resolve_fixed_ions(std::span<const ImportedMoleculeRecord> records,
                                      std::span<const std::string> component_ids)
    -> calculation::FixedIons;

} // namespace chargefw::io::detail
