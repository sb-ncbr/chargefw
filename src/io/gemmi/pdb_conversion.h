#pragma once

#include <chargefw/io/molecule_record.h>

#include <string>
#include <vector>

#include <gemmi/cifdoc.hpp>

namespace chargefw::io::gemmi::pdb_conversion {

struct ConvertedPdb {
    ::gemmi::cif::Document document;
    // The _atom_site.id of each calculated atom: atom_site_ids[conformer][atom].
    std::vector<std::vector<std::string>> atom_site_ids;
};

// Converts the PDB text a record was read from to mmCIF as Gemmi does: the whole structure, with
// its entities, secondary structure, assemblies, and LINK and SSBOND connections. The document
// shows what was calculated; it does not state the CONECT bonds or the bond orders of the record,
// so reading it does not give the record again.
[[nodiscard]] auto convert(const std::string& text, const ImportedMoleculeRecord& record)
    -> ConvertedPdb;

} // namespace chargefw::io::gemmi::pdb_conversion
