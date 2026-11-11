#include <chargefw/parameters/models/common_parameters.h>

#include "parameters/parameter_helpers.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace chargefw::parameters::detail {
namespace {

auto validate_parameter(const NamedParameter& parameter, const std::string_view context) -> void {
    const auto& [name, value] = parameter;

    if (name.empty()) {
        throw std::invalid_argument{std::string{context} + " parameter name must not be empty"};
    }

    if (!std::isfinite(value)) {
        throw std::invalid_argument{std::string{context} + " parameter '" + name +
                                    "' must be finite"};
    }
}

} // namespace

auto same_names(const std::span<const NamedParameter> first,
                const std::span<const NamedParameter> second) -> bool {
    return std::ranges::equal(first, second, [](const auto& lhs, const auto& rhs) -> bool {
        return lhs.name == rhs.name;
    });
}

auto validate_named_parameters(std::span<const NamedParameter> parameters,
                               const std::string_view context) -> void {
    for (auto first = parameters.begin(); first != parameters.end(); ++first) {
        validate_parameter(*first, context);

        const auto& first_name = first->name;

        for (auto second = std::next(first); second != parameters.end(); ++second) {
            const auto& second_name = second->name;

            if (first_name == second_name) {
                throw std::invalid_argument{"duplicate " + std::string{context} + " parameter '" +
                                            first_name + "'"};
            }
        }
    }
}

auto contains_named_parameter(std::span<const NamedParameter> parameters,
                              const std::string_view name) noexcept -> bool {
    return std::ranges::any_of(
        parameters, [name](const auto& parameter) -> bool { return parameter.name == name; });
}

auto named_parameter(std::span<const NamedParameter> parameters, const std::string_view name,
                     const std::string_view context) -> double {
    for (const auto& [parameter_name, value] : parameters) {
        if (parameter_name == name) {
            return value;
        }
    }

    throw std::invalid_argument{"unknown " + std::string{context} + " parameter '" +
                                std::string{name} + "'"};
}

} // namespace chargefw::parameters::detail
