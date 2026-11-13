#pragma once

#include "methods/method_applicability.h"

#include <chargefw/features/prepared_molecule_collection.h>

#include <optional>
#include <string>

namespace chargefw::methods::detail {

auto validate_selected_candidate(const ApplicableMethod& selected,
                                 const features::PreparedMoleculeCollection& molecules,
                                 bool require_fixed_ions = false) -> void;
auto validate_coordinate_targets(const ApplicableMethod& selected,
                                 const features::PreparedMoleculeCollection& molecules) -> void;
[[nodiscard]] auto parameter_set_id_for(const ApplicableMethod& selected)
    -> std::optional<std::string>;

} // namespace chargefw::methods::detail
