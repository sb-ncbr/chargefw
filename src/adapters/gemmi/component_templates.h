#pragma once

#include <chargefw/core/bond.h>

#include <cstdint>
#include <span>
#include <string_view>

namespace chargefw::adapters::gemmi::component_templates {

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

} // namespace chargefw::adapters::gemmi::component_templates
