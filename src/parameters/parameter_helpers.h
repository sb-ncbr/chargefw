#pragma once

#include <chargefw/parameters/models/common_parameters.h>

#include <span>
#include <string_view>

namespace chargefw::parameters::detail {

[[nodiscard]] auto same_names(std::span<const NamedParameter> first,
                              std::span<const NamedParameter> second) -> bool;

auto validate_named_parameters(std::span<const NamedParameter> parameters, std::string_view context)
    -> void;

[[nodiscard]] auto contains_named_parameter(std::span<const NamedParameter> parameters,
                                            std::string_view name) noexcept -> bool;

[[nodiscard]] auto named_parameter(std::span<const NamedParameter> parameters,
                                   std::string_view name, std::string_view context) -> double;

} // namespace chargefw::parameters::detail
