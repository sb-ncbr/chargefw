#pragma once

#include <chargefw/core/molecule.h>
#include <chargefw/core/periodic_table.h>

#include <algorithm>
#include <cstddef>
#include <map>
#include <span>
#include <string>
#include <utility>
#include <vector>

namespace chargefw::core::detail {

[[nodiscard]] inline auto bond_order_description(const BondOrder order) -> std::string_view {
    switch (order) {
    case BondOrder::SINGLE:
        return "single";
    case BondOrder::DOUBLE:
        return "double";
    case BondOrder::TRIPLE:
        return "triple";
    }
    return "unknown";
}

[[nodiscard]] inline auto atom_description(const Molecule& molecule, const std::size_t atom_index)
    -> std::string {
    const auto& atom = molecule.atom(atom_index);
    auto description = "atom " + std::to_string(atom_index + 1) + " (";
    if (!atom.name().empty()) {
        description += "source name '" + std::string{atom.name()} + "', ";
    }
    description += std::string{element_symbol(atom.atomic_number())} + ", formal charge " +
                   (atom.formal_charge() > 0 ? "+" : "") + std::to_string(atom.formal_charge()) +
                   ")";
    return description;
}

[[nodiscard]] inline auto bond_description(const Molecule& molecule, const std::size_t bond_index)
    -> std::string {
    const auto& bond = molecule.bond(bond_index);
    return "bond " + std::to_string(bond_index + 1) + " between " +
           atom_description(molecule, bond.first_atom_index()) + " and " +
           atom_description(molecule, bond.second_atom_index()) + ", " +
           std::string{bond_order_description(bond.order())} + " bond";
}

[[nodiscard]] inline auto molecule_description(const Molecule& molecule,
                                               const std::size_t molecule_index) -> std::string {
    auto description = "molecule " + std::to_string(molecule_index + 1);
    if (!molecule.name().empty()) {
        description += " ('" + std::string{molecule.name()} + "')";
    }
    return description;
}

// The atoms per element, most frequent first, such as "C: 2860, O: 797".
[[nodiscard]] inline auto element_counts_description(const Molecule& molecule,
                                                     const std::span<const std::size_t> atoms)
    -> std::string {
    auto counts = std::map<int, std::size_t>{};
    for (const auto atom_index : atoms) {
        ++counts[molecule.atom(atom_index).atomic_number()];
    }
    auto ordered = std::vector<std::pair<int, std::size_t>>{counts.begin(), counts.end()};
    std::ranges::stable_sort(ordered, std::ranges::greater{}, &std::pair<int, std::size_t>::second);
    auto description = std::string{};
    for (const auto& [atomic_number, count] : ordered) {
        if (!description.empty()) {
            description += ", ";
        }
        description += std::string{element_symbol(atomic_number)} + ": " + std::to_string(count);
    }
    return description;
}

} // namespace chargefw::core::detail
