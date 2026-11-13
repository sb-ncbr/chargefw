#pragma once

#include <chargefw/io/molecule_record.h>

#include <istream>
#include <optional>
#include <string>

namespace chargefw::io {

// Parses exactly one MOL record and throws on invalid input. The record ends at M  END; a following
// SDF separator is not consumed. V2000 atom/bond blocks plus M  CHG and the V3000 CTAB atom/bond
// subset are supported.
[[nodiscard]] auto parse_mol(std::istream& input, MoleculeRecordIdentity identity)
    -> ImportedMoleculeRecord;

// Stateful single-record counterpart to SdfReader. The first call parses the MOL record and later
// calls return clean end-of-file.
class MolReader {
  public:
    explicit MolReader(std::istream& input, std::string source_name = {});

    MolReader(const MolReader&) = delete;
    auto operator=(const MolReader&) -> MolReader& = delete;
    MolReader(MolReader&&) = default;
    auto operator=(MolReader&&) -> MolReader& = default;
    ~MolReader() = default;

    [[nodiscard]] auto next() -> std::optional<ImportedMoleculeRecord>;

  private:
    std::istream* input_;
    std::string source_;
    bool consumed_ = false;
};

} // namespace chargefw::io
