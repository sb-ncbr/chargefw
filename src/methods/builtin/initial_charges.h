#pragma once

#include <chargefw/core/molecule.h>
#include <chargefw/methods/method_options.h>

#include <cstddef>
#include <string>
#include <variant>
#include <vector>

namespace chargefw::methods::builtin::detail {

// Prerequisite checks also inspect options that have failed schema validation.
[[nodiscard]] inline auto uses_formal_initial_charges(const MethodOptions& options) -> bool {
    const auto option = options.values().find("initial_charges");
    if (option == options.values().end()) {
        return false;
    }
    const auto* value = std::get_if<std::string>(&option->second);
    return value != nullptr && *value == "formal";
}

[[nodiscard]] inline auto formal_initial_charges(const core::Molecule& molecule)
    -> std::vector<double> {
    auto charges = std::vector<double>(molecule.atom_count(), 0.0);
    for (std::size_t atom_index = 0; atom_index < molecule.atom_count(); ++atom_index) {
        charges[atom_index] = static_cast<double>(molecule.atom(atom_index).formal_charge());
    }
    return charges;
}

} // namespace chargefw::methods::builtin::detail
