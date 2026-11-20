#pragma once

#include <cstdint>
#include <stdexcept>
#include <string>
#include <string_view>

namespace chargefw::parameters {

// Why a parameter-set file was refused. More kinds may be added, so a switch needs a default
// branch.
enum class ParseErrorKind : std::uint8_t {
    // The file breaks the format, for example it is not JSON, lacks a member, or has a value of
    // the wrong type.
    malformed,
    // The file is in a format version this version of ChargeFW does not read.
    unsupported,
};

[[nodiscard]] auto to_string(ParseErrorKind value) -> std::string_view;

// A parameter-set file whose content the reader refuses. kind() gives the cause as a code and
// reason() as text, starting with the member of the document, such as
// "metadata.format_version: ...". source_name() names the file or stream, or is empty, and what()
// gives source name and reason. A file that cannot be opened is a std::runtime_error that is not a
// ParseError.
class ParseError : public std::runtime_error {
  public:
    ParseError(ParseErrorKind kind, std::string reason, std::string source_name = {});

    [[nodiscard]] auto kind() const noexcept -> ParseErrorKind;
    [[nodiscard]] auto reason() const noexcept -> const std::string&;
    [[nodiscard]] auto source_name() const noexcept -> const std::string&;

  private:
    ParseErrorKind kind_;
    std::string reason_;
    std::string source_name_;
};

} // namespace chargefw::parameters
