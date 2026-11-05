#pragma once

#include <chargefw/parameters/models/parameter_set.h>
#include <chargefw/parameters/models/parameter_set_metadata.h>

#include <cstddef>
#include <cstdint>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace chargefw::parameters {

enum class ParameterScope : std::uint8_t { common, atom, bond };

[[nodiscard]] auto to_string(ParameterScope value) -> std::string_view;

// One value of a parameter set: a common parameter, or one parameter of one atom or bond entry.
struct ParameterSlot {
    ParameterScope scope = ParameterScope::common;
    // Index of the atom or bond entry; 0 for common parameters.
    std::size_t entry_index = 0;
    std::string name;

    [[nodiscard]] auto operator==(const ParameterSlot&) const -> bool = default;
};

// Every value of `set` in a fixed order: the common parameters, then each atom entry with its
// parameters, then each bond entry with its parameters, each in stored order.
[[nodiscard]] auto parameter_layout(const ParameterSet& set) -> std::vector<ParameterSlot>;

// The values of `set` in the order of parameter_layout().
[[nodiscard]] auto parameter_values(const ParameterSet& set) -> std::vector<double>;

// True when `candidate` has the method ID, entry keys, and parameter names of `reference`, in the
// same order; the values and the descriptive metadata may differ.
[[nodiscard]] auto compatible(const ParameterSet& reference, const ParameterSet& candidate) -> bool;

// A copy of `set` with `values` in the order of parameter_layout() and the identifying `metadata`.
// The method ID is taken from `set`; metadata.method_id must be empty or equal to it. Throws
// std::invalid_argument for a wrong number of values, a non-finite value, or another method ID.
[[nodiscard]] auto with_parameter_values(const ParameterSet& set, std::span<const double> values,
                                         ParameterSetMetadata metadata) -> ParameterSet;

} // namespace chargefw::parameters
