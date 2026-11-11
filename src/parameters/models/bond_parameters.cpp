#include <chargefw/parameters/models/bond_parameters.h>

#include "parameters/parameter_helpers.h"

#include <utility>

namespace chargefw::parameters {
namespace {

auto validate_entries(std::span<const BondParameterEntry> entries) -> void {
    for (const auto& [key, parameters] : entries) {
        detail::validate_named_parameters(parameters, "bond");
    }
}

} // namespace

BondParameters::BondParameters(std::vector<BondParameterEntry> entries)
    : entries_{std::move(entries)} {
    validate_entries(entries_);
}

auto BondParameters::entries() const noexcept -> std::span<const BondParameterEntry> {
    return entries_;
}

auto BondParameters::size() const noexcept -> std::size_t {
    return entries_.size();
}

auto BondParameters::empty() const noexcept -> bool {
    return entries_.empty();
}

auto BondParameters::operator[](const std::size_t index) const noexcept
    -> const BondParameterEntry& {
    return entries_[index];
}

auto BondParameters::at(const std::size_t index) const -> const BondParameterEntry& {
    return entries_.at(index);
}

auto BondParameters::contains(const std::size_t entry_index, const std::string_view name) const
    -> bool {
    return detail::contains_named_parameter(at(entry_index).parameters, name);
}

auto BondParameters::parameter(const std::size_t entry_index, const std::string_view name) const
    -> double {
    return detail::named_parameter(at(entry_index).parameters, name, "bond");
}

} // namespace chargefw::parameters
