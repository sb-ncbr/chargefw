#pragma once

#include <chargefw/adapters/gemmi/input_options.h>
#include <chargefw/adapters/molecule_record.h>

#include <istream>
#include <optional>
#include <string>

namespace chargefw::adapters::gemmi::pdb_input {

// Reads a PDB structure through Gemmi. The first PDB MODEL defines atom topology and each model
// becomes a conformer after validating the same atom sequence. Record selection, alternate-location
// selection (the first location in source order), and bond strategies follow the PDB and mmCIF
// input contract in docs/reference/formats.md.
class PdbReader {
  public:
    explicit PdbReader(std::istream& input, std::string source = {},
                       ::chargefw::adapters::gemmi::InputOptions options = {});

    PdbReader(const PdbReader&) = delete;
    auto operator=(const PdbReader&) -> PdbReader& = delete;
    PdbReader(PdbReader&&) = default;
    auto operator=(PdbReader&&) -> PdbReader& = default;
    ~PdbReader() = default;

    [[nodiscard]] auto next() -> std::optional<ImportedMoleculeRecord>;
    [[nodiscard]] auto options() const noexcept -> ::chargefw::adapters::gemmi::InputOptions;

  private:
    std::optional<ImportedMoleculeRecord> record_;
    ::chargefw::adapters::gemmi::InputOptions options_;
};

} // namespace chargefw::adapters::gemmi::pdb_input
