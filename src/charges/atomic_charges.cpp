#include <chargefw/charges/atomic_charges.h>

#include <cmath>
#include <numeric>
#include <stdexcept>
#include <utility>

namespace chargefw::charges {
namespace {

auto validate_charges(std::span<const double> values) -> void {
    for (const auto value : values) {
        if (!std::isfinite(value)) {
            throw NonFiniteCharges{"atomic charges must be finite"};
        }
    }
}

} // namespace

AtomicCharges::AtomicCharges(std::vector<double> values) : values_{std::move(values)} {
    validate_charges(values_);
}

auto AtomicCharges::values() const noexcept -> std::span<const double> {
    return values_;
}

auto AtomicCharges::size() const noexcept -> std::size_t {
    return values_.size();
}

auto AtomicCharges::empty() const noexcept -> bool {
    return values_.empty();
}

auto AtomicCharges::total() const noexcept -> double {
    return std::accumulate(values_.begin(), values_.end(), 0.0);
}

[[nodiscard]] auto AtomicCharges::operator[](std::size_t index) const noexcept -> double {
    return values_[index];
}

[[nodiscard]] auto AtomicCharges::at(std::size_t index) const -> double {
    return values_.at(index);
}

} // namespace chargefw::charges