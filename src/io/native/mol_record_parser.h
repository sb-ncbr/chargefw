#pragma once

#include <chargefw/io/molecule_record.h>

#include <cstddef>
#include <istream>

namespace chargefw::io::native::mol_input {

// Internal overload for an SDF reader's running source-line counter.
auto parse_mol(std::istream& input, MoleculeRecordIdentity identity, std::size_t& line)
    -> ImportedMoleculeRecord;

} // namespace chargefw::io::native::mol_input
