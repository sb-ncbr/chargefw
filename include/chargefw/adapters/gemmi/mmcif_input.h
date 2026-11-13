#pragma once

#include <chargefw/adapters/gemmi/input_options.h>
#include <chargefw/adapters/molecule_record.h>

#include <istream>
#include <memory>
#include <optional>
#include <string>

namespace chargefw::adapters::gemmi::mmcif_input {

// Parses an mmCIF document eagerly through Gemmi, then converts one coordinate-bearing data block
// per next() call. Each such block becomes one record; models within a block become conformers
// after atom-sequence validation. Selection, alternate-location, and bond-strategy behavior match
// pdb_input::PdbReader.
class MmcifReader {
  public:
    explicit MmcifReader(std::istream& input, std::string source = {},
                         ::chargefw::adapters::gemmi::InputOptions options = {});

    MmcifReader(const MmcifReader&) = delete;
    auto operator=(const MmcifReader&) -> MmcifReader& = delete;
    MmcifReader(MmcifReader&&) noexcept;
    auto operator=(MmcifReader&&) noexcept -> MmcifReader&;
    ~MmcifReader();

    [[nodiscard]] auto next() -> std::optional<ImportedMoleculeRecord>;

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace chargefw::adapters::gemmi::mmcif_input
