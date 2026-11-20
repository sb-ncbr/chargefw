#pragma once

#include <chargefw/core/issue_location.h>

#include <cstddef>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>

namespace chargefw::io {

// Why a writer refused a result. More kinds may be added, so a switch needs a default branch.
enum class OutputErrorKind : std::uint8_t {
    // The format needs coordinates and a molecule has none.
    missing_coordinates,
    // The format cannot hold the number of results, molecules, conformers, or atoms, such as
    // several molecules in PQR or a molecule without atoms.
    unsupported_shape,
    // The format cannot hold a value, such as a chain ID of two characters in PQR.
    unrepresentable,
    // The document or toolkit molecule to add charges to is not the one that was read.
    target_mismatch,
    // The target already has charges that were not to be replaced, or that cannot be added to.
    existing_charges,
    // The calculation failed, so there are no charges to write.
    no_charges,
    // The format is written only for molecules read from other input, such as PQR, which needs a
    // structure read from PDB or mmCIF.
    unsupported_source,
};

[[nodiscard]] auto to_string(OutputErrorKind value) -> std::string_view;

// A result that a writer cannot write in its format. kind() gives the cause as a code and reason()
// as text, with the value and what the format requires. molecule_index() is the input molecule the
// refusal concerns, counted from zero, or nothing when it concerns the result as a whole, such as
// too many molecules for the format; location() gives the conformer, atoms, and bonds. what() gives
// where and why, counted from one, such as "molecule 1, atom 5: PQR output cannot write chain ID
// 'AB': ...".
class OutputError : public std::invalid_argument {
  public:
    OutputError(OutputErrorKind kind, std::string reason,
                std::optional<std::size_t> molecule_index = {}, core::IssueLocation location = {});

    [[nodiscard]] auto kind() const noexcept -> OutputErrorKind;
    [[nodiscard]] auto reason() const noexcept -> const std::string&;
    [[nodiscard]] auto molecule_index() const noexcept -> std::optional<std::size_t>;
    [[nodiscard]] auto location() const noexcept -> const core::IssueLocation&;

  private:
    OutputErrorKind kind_;
    std::string reason_;
    std::optional<std::size_t> molecule_index_;
    core::IssueLocation location_;
};

} // namespace chargefw::io
