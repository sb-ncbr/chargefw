#pragma once

#include <chargefw/core/bond.h>

#include <cstdint>
#include <functional>
#include <optional>
#include <span>
#include <string_view>

namespace chargefw::io::gemmi::component_templates {

struct BondTemplate {
    std::string_view first;
    std::string_view second;
    core::BondOrder order;
};

enum class ComponentKind : std::uint8_t {
    amino_acid,
    nucleotide,
    water,
};

// A group whose Lewis structure depends on where the residue is protonated, such as the histidine
// ring. Each alternative lists the bonds of the group with orders other than the template's.
struct ProtonationMotif {
    std::span<const std::string_view> atoms;
    std::span<const std::span<const BondTemplate>> alternatives;
};

struct ComponentTemplate {
    ComponentKind kind;
    std::span<const BondTemplate> bonds;
    const ProtonationMotif* motif = nullptr;
};

[[nodiscard]] auto find(std::string_view component) -> const ComponentTemplate*;

struct MotifAtom {
    int atomic_number;
    int formal_charge;
};

// The atom of a residue with the given name, or nothing when the residue has none.
using MotifAtomLookup = std::function<std::optional<MotifAtom>(std::string_view)>;

// The bonds whose orders replace the template's for one residue: those of the first alternative
// Lewis structure in which every motif atom has the bond-order sum its element and formal charge
// require, unless the template already has it. A residue that no structure fits, such as one
// missing a motif atom, keeps the template orders.
[[nodiscard]] auto protonation_bonds(const ComponentTemplate& component_template,
                                     const MotifAtomLookup& atom) -> std::span<const BondTemplate>;

// Whether the order of a template bond depends on the residue's protonation state.
[[nodiscard]] auto varies_with_protonation(const ComponentTemplate& component_template,
                                           std::string_view first, std::string_view second) -> bool;

// The order of a template bond in a residue with the given protonation bonds, or nothing when the
// template has no such bond.
[[nodiscard]] auto bond_order(const ComponentTemplate& component_template,
                              std::span<const BondTemplate> protonation, std::string_view first,
                              std::string_view second) -> std::optional<core::BondOrder>;

} // namespace chargefw::io::gemmi::component_templates
