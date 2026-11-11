#include <chargefw/parameters/models/atom_parameters.h>

#include "parameters/parameter_helpers.h"

#include <utility>

namespace chargefw::parameters {
namespace {

auto validate_entries(std::span<const AtomParameterEntry> entries) -> void {
    for (const auto& [key, parameters] : entries) {
        detail::validate_named_parameters(parameters, "atom");
    }
}

} // namespace

AtomParameters::AtomParameters(std::vector<AtomParameterEntry> entries)
    : entries_{std::move(entries)} {
    validate_entries(entries_);
}

auto AtomParameters::entries() const noexcept -> std::span<const AtomParameterEntry> {
    return entries_;
}

auto AtomParameters::size() const noexcept -> std::size_t {
    return entries_.size();
}

auto AtomParameters::empty() const noexcept -> bool {
    return entries_.empty();
}

auto AtomParameters::operator[](const std::size_t index) const noexcept
    -> const AtomParameterEntry& {
    return entries_[index];
}

auto AtomParameters::at(const std::size_t index) const -> const AtomParameterEntry& {
    return entries_.at(index);
}

auto AtomParameters::contains(const std::size_t entry_index, const std::string_view name) const
    -> bool {
    return detail::contains_named_parameter(at(entry_index).parameters, name);
}

auto AtomParameters::parameter(const std::size_t entry_index, const std::string_view name) const
    -> double {
    return detail::named_parameter(at(entry_index).parameters, name, "atom");
}

} // namespace chargefw::parameters
