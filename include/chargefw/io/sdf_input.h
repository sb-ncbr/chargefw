#pragma once

#include <chargefw/io/molecule_record.h>

#include <cstddef>
#include <istream>
#include <optional>
#include <string>

namespace chargefw::io {

// Bounded-memory SDF reader. Each call consumes at most one SDF record. A null result denotes clean
// end-of-file after at least one record; input with no record, and invalid input, throw ParseError.
class SdfReader {
  public:
    explicit SdfReader(std::istream& input, std::string source_name = {});

    SdfReader(const SdfReader&) = delete;
    auto operator=(const SdfReader&) -> SdfReader& = delete;
    SdfReader(SdfReader&&) = default;
    auto operator=(SdfReader&&) -> SdfReader& = default;
    ~SdfReader() = default;

    [[nodiscard]] auto next() -> std::optional<ImportedMoleculeRecord>;

  private:
    std::istream* input_;
    std::string source_;
    std::size_t record_index_ = 0;
    std::size_t line_ = 0;
};

} // namespace chargefw::io
