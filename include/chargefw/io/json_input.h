#pragma once

#include <chargefw/io/conformer_selection.h>
#include <chargefw/io/molecule_record.h>

#include <istream>
#include <memory>
#include <optional>
#include <string>

namespace chargefw::io {

// ChargeFW JSON reader for a version 1.0 molecule document. Each call returns the next molecule
// record from the top-level `molecules` array. A null result denotes clean end-of-document; invalid
// input throws and terminates import. Molecule array order is source atom/conformer order; atom
// names are intentionally unsupported.
class JsonReader {
  public:
    explicit JsonReader(std::istream& input, std::string source_name = {},
                        ConformerSelection conformers = ConformerSelection::all);

    JsonReader(const JsonReader&) = delete;
    auto operator=(const JsonReader&) -> JsonReader& = delete;
    JsonReader(JsonReader&&) noexcept;
    auto operator=(JsonReader&&) noexcept -> JsonReader&;
    ~JsonReader();

    [[nodiscard]] auto next() -> std::optional<ImportedMoleculeRecord>;

  private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace chargefw::io
