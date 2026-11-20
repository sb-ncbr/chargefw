#pragma once

#include <chargefw/io/molecule_record.h>

#include <cstddef>
#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>

namespace chargefw::io {

// Where in the source a reader stopped, as far as it knows.
struct ParseLocation {
    std::string source_name;
    std::optional<std::size_t> line = {};         // One-based line in the source.
    std::optional<std::size_t> record_index = {}; // Zero-based record in the source.
    PortableId record_id = {};
};

// Why a reader refused its input. More kinds may be added, so a switch needs a default branch.
enum class ParseErrorKind : std::uint8_t {
    // Nothing to read: an empty file, or a file of another format.
    no_records,
    // Text the reader cannot decode, such as a character that is not ASCII.
    invalid_encoding,
    // The file breaks its format, such as wrong counts or a bond to an atom that does not exist.
    malformed,
    // Valid data that ChargeFW does not read, such as query atoms.
    unsupported,
    // The models of a PDB or mmCIF structure do not have the same atoms.
    inconsistent_models,
};

[[nodiscard]] auto to_string(ParseErrorKind value) -> std::string_view;

// Input that a reader cannot read. kind() gives the cause as a code and reason() as text;
// location() gives where, and what() gives reason and location, such as
// "input.sdf:24: record 3: unsupported MOL version".
class ParseError : public std::runtime_error {
  public:
    ParseError(ParseErrorKind kind, std::string reason, ParseLocation location);

    [[nodiscard]] auto kind() const noexcept -> ParseErrorKind;
    [[nodiscard]] auto reason() const noexcept -> const std::string&;
    [[nodiscard]] auto location() const noexcept -> const ParseLocation&;

  private:
    ParseErrorKind kind_;
    std::string reason_;
    ParseLocation location_;
};

} // namespace chargefw::io
