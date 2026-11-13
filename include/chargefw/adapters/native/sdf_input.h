#pragma once

#include <chargefw/adapters/molecule_record.h>

#include <cstddef>
#include <istream>
#include <optional>
#include <string>

namespace chargefw::adapters {

// Bounded-memory SDF reader. Each call consumes at most one SDF record. A null result denotes clean
// end-of-file; invalid input throws and terminates import.
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

} // namespace chargefw::adapters
