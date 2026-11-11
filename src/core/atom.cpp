#include <chargefw/core/atom.h>
#include <chargefw/core/periodic_table.h>

#include <stdexcept>
#include <string>
#include <utility>

namespace chargefw::core {
namespace {

auto validate_atomic_number(int atomic_number) -> void {
    if (!periodic_table().contains(atomic_number)) {
        throw std::invalid_argument{"atomic number must be in range 1.." +
                                    std::to_string(periodic_table().elements().size())};
    }
}

} // namespace

Atom::Atom(int atomic_number, int formal_charge, std::string source_name)
    : atomic_number_{atomic_number}, formal_charge_{formal_charge}, name_{std::move(source_name)} {
    validate_atomic_number(atomic_number_);
}

auto Atom::atomic_number() const noexcept -> int {
    return atomic_number_;
}

auto Atom::formal_charge() const noexcept -> int {
    return formal_charge_;
}

auto Atom::name() const noexcept -> std::string_view {
    return name_;
}

} // namespace chargefw::core
