#pragma once

#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/conformer.h>
#include <chargefw/io/molecule_record.h>

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace chargefw::io {

// What every reader builds its records with.

// The warning that a record's aromatic bonds, given by their indices in the record in bond order,
// were read as single bonds.
[[nodiscard]] auto aromatic_bonds_issue(std::vector<std::size_t> bonds) -> ImportIssue;

// A record named after its record ID when `name` is empty.
[[nodiscard]] auto make_record(std::vector<core::Atom> atoms, std::vector<core::Bond> bonds,
                               std::vector<core::Conformer> conformers,
                               MoleculeRecordIdentity identity, std::string name = {},
                               std::vector<ImportIssue> issues = {},
                               std::optional<MoleculeImportMetadata> import_metadata = std::nullopt)
    -> ImportedMoleculeRecord;

} // namespace chargefw::io
