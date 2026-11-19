#pragma once

#include <chargefw/core/issue_location.h>

#include <nanobind/nanobind.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/string.h>
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

// An issue with its location, as the Python issue types read it; every issue type has the same
// fields. The kind's to_string() is found in the namespace of the kind.
template <typename Issue> auto issue_payload(const Issue& issue) -> nanobind::dict {
    auto result = nanobind::dict{};
    result["kind"] = std::string{to_string(issue.kind)};
    result["message"] = issue.message;
    result["molecule_index"] = issue.molecule_index;
    add_issue_location(result, issue.location);
    return result;
}

} // namespace chargefw::python
