#pragma once

#include <chargefw/core/issue_location.h>

#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>

namespace chargefw::calculation {

// A calculated charge is probably wrong when it is further than this from the atom's formal
// charge. Recommended methods on clean inputs stay within 4.9 e, the charge of QEq phosphorus in
// DNA; broken inputs and method breakdowns reach 6.1 e or more.
inline constexpr double charge_deviation_limit = 6.0; // e

enum class ChargeIssueKind : std::uint8_t {
    // Charges are further than charge_deviation_limit from their formal charges.
    charge_out_of_range,
};

[[nodiscard]] auto to_string(ChargeIssueKind value) -> std::string_view;

// A caution about calculated charges: they exist, but are probably wrong. Execution finds them
// for the charges of the executed plan, one issue per assignment at most.
struct ChargeIssue {
    ChargeIssueKind kind;
    std::string message;
    std::size_t molecule_index; // Zero-based in the supplied molecule collection.
    // The conformer of the charges, empty for methods that do not use coordinates; the atom the
    // message names, the one furthest from its formal charge; and every affected atom.
    core::IssueLocation location = {};
};

} // namespace chargefw::calculation
