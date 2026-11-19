#pragma once

#include <chargefw/core/issue_location.h>

#include <nanobind/nanobind.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/vector.h>

namespace chargefw::python {

// The location fields that every Python issue type has.
inline auto add_issue_location(nanobind::dict& payload, const core::IssueLocation& location)
    -> void {
    payload["conformer_index"] = location.conformer_index;
    payload["atom_index"] = location.atom_index;
    payload["other_atom_index"] = location.other_atom_index;
    payload["bond_index"] = location.bond_index;
    payload["atom_indices"] = location.atom_indices;
    payload["bond_indices"] = location.bond_indices;
}

} // namespace chargefw::python
