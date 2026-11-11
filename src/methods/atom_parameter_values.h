#pragma once

#include <chargefw/parameters/models/parameter_view.h>

#include <cstddef>
#include <string_view>
#include <vector>

namespace chargefw::methods::detail {

[[nodiscard]] inline auto atom_parameter_values(const parameters::ParameterView& view,
                                                const std::string_view name,
                                                const std::size_t atom_count)
    -> std::vector<double> {
    const auto parameter = view.atom(name);
    auto values = std::vector<double>(atom_count);
    for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
        values[atom_index] = parameter[atom_index];
    }
    return values;
}

} // namespace chargefw::methods::detail
