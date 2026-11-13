#pragma once

#include <chargefw/core/molecule.h>

#include <string_view>
#include <vector>

namespace chargefw::io::native::detail {

// SYBYL atom types for MOL2 output, in atom order, derived from elements, bonds, bond orders, and
// formal charges. Aromaticity is not perceived: ring atoms of a Kekulé structure get their sp2
// types, such as C.2. Elements without a hybridization-specific type get their element symbol.
[[nodiscard]] auto sybyl_atom_types(const core::Molecule& molecule)
    -> std::vector<std::string_view>;

} // namespace chargefw::io::native::detail
