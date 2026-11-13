#include <chargefw/charges/charge_set.h>

#include <stdexcept>
#include <utility>

namespace chargefw::charges {
namespace {

auto validate_method_id(std::string_view method_id) -> void {
    if (method_id.empty()) {
        throw std::invalid_argument{"method id must not be empty"};
    }
}

} // namespace

ChargeSet::ChargeSet(std::string method_id, std::vector<ChargeAssignment> assignments,
                     std::optional<std::string> parameter_set_id)
    : method_id_{std::move(method_id)}, parameter_set_id_{std::move(parameter_set_id)},
      assignments_{std::move(assignments)} {
    validate_method_id(method_id_);
}

auto ChargeSet::method_id() const noexcept -> std::string_view {
    return method_id_;
}

auto ChargeSet::parameter_set_id() const noexcept -> std::optional<std::string_view> {
    if (!parameter_set_id_.has_value()) {
        return std::nullopt;
    }

    return std::string_view{*parameter_set_id_};
}

auto ChargeSet::assignments() const noexcept -> std::span<const ChargeAssignment> {
    return assignments_;
}

auto ChargeSet::assignment(const std::size_t index) const -> const ChargeAssignment& {
    return assignments_.at(index);
}

auto ChargeSet::size() const noexcept -> std::size_t {
    return assignments_.size();
}

auto ChargeSet::empty() const noexcept -> bool {
    return assignments_.empty();
}

} // namespace chargefw::charges
