#include <chargefw/parameterization/parameter_vector.h>

#include <chargefw/parameters/models/atom_parameters.h>
#include <chargefw/parameters/models/bond_parameters.h>
#include <chargefw/parameters/models/common_parameters.h>

#include <cstddef>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace chargefw::parameterization {
namespace {

using parameters::NamedParameter;
using parameters::ParameterSet;
using parameters::ParameterSetMetadata;

[[nodiscard]] auto value_count(const ParameterSet& set) -> std::size_t {
    auto count = set.common().size();
    for (const auto& entry : set.atom().entries()) {
        count += entry.parameters.size();
    }
    for (const auto& entry : set.bond().entries()) {
        count += entry.parameters.size();
    }
    return count;
}

template <typename Entry>
auto append_slots(const std::span<const Entry> entries, const ParameterScope scope,
                  std::vector<ParameterSlot>& slots) -> void {
    for (std::size_t entry_index = 0; entry_index < entries.size(); ++entry_index) {
        for (const auto& parameter : entries[entry_index].parameters) {
            slots.push_back({.scope = scope, .entry_index = entry_index, .name = parameter.name});
        }
    }
}

// Copies `parameters` with their values replaced by the next values of `values`.
[[nodiscard]] auto take_values(const std::span<const NamedParameter> parameters,
                               std::span<const double>& values) -> std::vector<NamedParameter> {
    auto result = std::vector<NamedParameter>{};
    result.reserve(parameters.size());
    for (const auto& parameter : parameters) {
        result.push_back({.name = parameter.name, .value = values.front()});
        values = values.subspan(1);
    }
    return result;
}

template <typename Entry>
[[nodiscard]] auto take_entry_values(const std::span<const Entry> entries,
                                     std::span<const double>& values) -> std::vector<Entry> {
    auto result = std::vector<Entry>{};
    result.reserve(entries.size());
    for (const auto& entry : entries) {
        result.push_back({.key = entry.key, .parameters = take_values(entry.parameters, values)});
    }
    return result;
}

} // namespace

auto to_string(const ParameterScope value) -> std::string_view {
    switch (value) {
    case ParameterScope::common:
        return "common";
    case ParameterScope::atom:
        return "atom";
    case ParameterScope::bond:
        return "bond";
    }
    throw std::invalid_argument{"unknown parameter scope"};
}

auto parameter_layout(const ParameterSet& set) -> std::vector<ParameterSlot> {
    auto slots = std::vector<ParameterSlot>{};
    slots.reserve(value_count(set));
    for (const auto& parameter : set.common().parameters()) {
        slots.push_back(
            {.scope = ParameterScope::common, .entry_index = 0, .name = parameter.name});
    }
    append_slots(set.atom().entries(), ParameterScope::atom, slots);
    append_slots(set.bond().entries(), ParameterScope::bond, slots);
    return slots;
}

auto parameter_values(const ParameterSet& set) -> std::vector<double> {
    auto values = std::vector<double>{};
    values.reserve(value_count(set));
    for (const auto& parameter : set.common().parameters()) {
        values.push_back(parameter.value);
    }
    for (const auto& entry : set.atom().entries()) {
        for (const auto& parameter : entry.parameters) {
            values.push_back(parameter.value);
        }
    }
    for (const auto& entry : set.bond().entries()) {
        for (const auto& parameter : entry.parameters) {
            values.push_back(parameter.value);
        }
    }
    return values;
}

auto with_parameter_values(const ParameterSet& set, std::span<const double> values,
                           ParameterSetMetadata metadata) -> ParameterSet {
    const auto expected = value_count(set);
    if (values.size() != expected) {
        throw std::invalid_argument{"parameter set '" + std::string{set.id()} + "' has " +
                                    std::to_string(expected) + " values, but " +
                                    std::to_string(values.size()) + " were given"};
    }
    if (!metadata.method_id.empty() && metadata.method_id != set.method_id()) {
        throw std::invalid_argument{"parameter values of method '" + std::string{set.method_id()} +
                                    "' cannot form a set of method '" + metadata.method_id + "'"};
    }
    metadata.method_id = set.method_id();

    // The constructors reject non-finite values.
    auto common = parameters::CommonParameters{take_values(set.common().parameters(), values)};
    auto atom = parameters::AtomParameters{take_entry_values(set.atom().entries(), values)};
    auto bond = parameters::BondParameters{take_entry_values(set.bond().entries(), values)};
    return ParameterSet{std::move(metadata), std::move(common), std::move(atom), std::move(bond)};
}

} // namespace chargefw::parameterization
