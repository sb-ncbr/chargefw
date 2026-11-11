#include <chargefw/parameters/models/parameter_set.h>

#include "parameters/parameter_helpers.h"

#include <algorithm>
#include <span>
#include <stdexcept>
#include <utility>

namespace chargefw::parameters {
namespace {

template <typename Entry>
[[nodiscard]] auto same_entries(const std::span<const Entry> first,
                                const std::span<const Entry> second) -> bool {
    return std::ranges::equal(first, second, [](const Entry& lhs, const Entry& rhs) -> bool {
        return lhs.key == rhs.key && detail::same_names(lhs.parameters, rhs.parameters);
    });
}

auto validate_metadata(const ParameterSetMetadata& metadata) -> void {
    if (metadata.id.empty()) {
        throw std::invalid_argument{"parameter set id must not be empty"};
    }

    if (metadata.method_id.empty()) {
        throw std::invalid_argument{"parameter set method id must not be empty"};
    }
}

} // namespace

ParameterSet::ParameterSet(ParameterSetMetadata metadata, CommonParameters common,
                           AtomParameters atom, BondParameters bond)
    : metadata_{std::move(metadata)}, common_{std::move(common)}, atom_{std::move(atom)},
      bond_{std::move(bond)} {
    validate_metadata(metadata_);
}

auto ParameterSet::metadata() const noexcept -> const ParameterSetMetadata& {
    return metadata_;
}

auto ParameterSet::id() const noexcept -> std::string_view {
    return metadata_.id;
}

auto ParameterSet::method_id() const noexcept -> std::string_view {
    return metadata_.method_id;
}

auto ParameterSet::name() const noexcept -> std::string_view {
    return metadata_.name;
}

auto ParameterSet::priority() const noexcept -> std::uint16_t {
    return metadata_.priority;
}

auto ParameterSet::common() const noexcept -> const CommonParameters& {
    return common_;
}

auto ParameterSet::atom() const noexcept -> const AtomParameters& {
    return atom_;
}

auto ParameterSet::bond() const noexcept -> const BondParameters& {
    return bond_;
}

auto compatible(const ParameterSet& reference, const ParameterSet& candidate) -> bool {
    return reference.method_id() == candidate.method_id() &&
           detail::same_names(reference.common().parameters(), candidate.common().parameters()) &&
           same_entries(reference.atom().entries(), candidate.atom().entries()) &&
           same_entries(reference.bond().entries(), candidate.bond().entries());
}

} // namespace chargefw::parameters
