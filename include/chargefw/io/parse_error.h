#pragma once

#include <chargefw/io/molecule_record.h>

#include <cstddef>
#include <optional>
#include <stdexcept>
#include <string>

namespace chargefw::io {

// Where in the source a reader stopped, as far as it knows.
struct ParseLocation {
    std::string source_name;
    std::optional<std::size_t> line = {};         // One-based line in the source.
    std::optional<std::size_t> record_index = {}; // Zero-based record in the source.
    PortableId record_id = {};
};

// Input that a reader cannot read. reason() gives only why; location() gives where, and what()
// gives both, such as "input.sdf:24: record 3: unsupported MOL version".
class ParseError : public std::runtime_error {
  public:
    ParseError(std::string reason, ParseLocation location);

    [[nodiscard]] auto reason() const noexcept -> const std::string&;
    [[nodiscard]] auto location() const noexcept -> const ParseLocation&;

  private:
    std::string reason_;
    ParseLocation location_;
};

} // namespace chargefw::io
