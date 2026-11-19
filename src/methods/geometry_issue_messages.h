#pragma once

#include "core/diagnostic_description.h"

#include <chargefw/core/issue_location.h>
#include <chargefw/core/molecule.h>

#include <string>

namespace chargefw::methods::detail {

// Why the coordinates of a conformer are unusable: two coincident atoms when the location has the
// other atom, otherwise one atom with non-finite coordinates. The conformer is in the location, not
// in the message. Fixed-ion planning rebuilds it for the original molecule, so it depends only on
// the molecule and the location.
[[nodiscard]] inline auto invalid_geometry_message(const core::Molecule& molecule,
                                                   const core::IssueLocation& location)
    -> std::string {
    if (location.other_atom_index.has_value()) {
        return core::detail::atom_description(molecule, *location.other_atom_index) + " and " +
               core::detail::atom_description(molecule, *location.atom_index) +
               " have identical coordinates; check their positions in the input structure";
    }
    return core::detail::atom_description(molecule, *location.atom_index) +
           " has non-finite coordinates; check its position in the input structure";
}

} // namespace chargefw::methods::detail
