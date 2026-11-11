#pragma once

#include <chargefw/adapters/gemmi/input_options.h>
#include <chargefw/adapters/molecule_record.h>

#include <gemmi/cifdoc.hpp>

#include <cstddef>
#include <istream>
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
    MmcifReader(MmcifReader&&) = default;
    auto operator=(MmcifReader&&) -> MmcifReader& = default;
    ~MmcifReader() = default;

    [[nodiscard]] auto next() -> std::optional<ImportedMoleculeRecord>;

  private:
    ::gemmi::cif::Document document_;
    std::string source_;
    ::chargefw::adapters::gemmi::InputOptions options_;
    std::size_t block_index_ = 0;
    std::size_t record_index_ = 0;
};

} // namespace chargefw::adapters::gemmi::mmcif_input
