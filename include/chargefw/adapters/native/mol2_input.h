#pragma once

#include <chargefw/adapters/molecule_record.h>

#include <cstddef>
#include <istream>
#include <optional>
#include <string>

namespace chargefw::adapters {

// Bounded-memory Tripos MOL2 reader. Each call consumes one @<TRIPOS>MOLECULE record. The reader
// uses the ATOM and BOND sections, with standard atom types and numeric bonds, and formal charges
// from UNITY_ATOM_ATTR; it skips other sections. Aromatic bond types are imported as single bonds.
class Mol2Reader {
  public:
    explicit Mol2Reader(std::istream& input, std::string source_name = {});

    Mol2Reader(const Mol2Reader&) = delete;
    auto operator=(const Mol2Reader&) -> Mol2Reader& = delete;
    Mol2Reader(Mol2Reader&&) = default;
    auto operator=(Mol2Reader&&) -> Mol2Reader& = default;
    ~Mol2Reader() = default;

    [[nodiscard]] auto next() -> std::optional<ImportedMoleculeRecord>;

  private:
    std::istream* input_;
    std::string source_;
    std::size_t record_index_ = 0;
    std::size_t line_ = 0;
    // The previous record ended at the MOLECULE line of the next one.
    bool next_record_started_ = false;
};

} // namespace chargefw::adapters
