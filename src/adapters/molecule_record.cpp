#include <chargefw/adapters/molecule_record.h>

#include "imported_records.h"

#include <chargefw/core/molecule.h>

#include <initializer_list>
#include <stdexcept>
#include <string>
#include <string_view>

namespace chargefw::adapters {

auto to_string(const MolecularSourceFormat value) -> std::string_view {
    switch (value) {
    case MolecularSourceFormat::molecule_json:
        return "molecule-json";
    case MolecularSourceFormat::mol:
        return "mol";
    case MolecularSourceFormat::sdf:
        return "sdf";
    case MolecularSourceFormat::mol2:
        return "mol2";
    case MolecularSourceFormat::pdb:
        return "pdb";
    case MolecularSourceFormat::mmcif:
        return "mmcif";
    }
    throw std::invalid_argument{"unknown molecular source format"};
}

auto molecular_source_format_from_string(const std::string_view value) -> MolecularSourceFormat {
    for (const auto format : {MolecularSourceFormat::molecule_json, MolecularSourceFormat::mol,
                              MolecularSourceFormat::sdf, MolecularSourceFormat::mol2,
                              MolecularSourceFormat::pdb, MolecularSourceFormat::mmcif}) {
        if (to_string(format) == value) {
            return format;
        }
    }
    throw std::invalid_argument{"unknown molecular source format: " + std::string{value}};
}

auto to_string(const SourceConnectivity value) -> std::string_view {
    switch (value) {
    case SourceConnectivity::absent:
        return "absent";
    case SourceConnectivity::explicitly_empty:
        return "explicitly-empty";
    case SourceConnectivity::present:
        return "present";
    }
    throw std::invalid_argument{"unknown source connectivity state"};
}

auto to_string(const ImportIssueKind value) -> std::string_view {
    switch (value) {
    case ImportIssueKind::partial_charges_ignored:
        return "partial_charges_ignored";
    case ImportIssueKind::coordinates_2d_ignored:
        return "coordinates_2d_ignored";
    case ImportIssueKind::v3000_bond_order_ignored:
        return "v3000_bond_order_ignored";
    case ImportIssueKind::aromatic_bonds_converted:
        return "aromatic_bonds_converted";
    case ImportIssueKind::symmetry_connections_ignored:
        return "symmetry_connections_ignored";
    case ImportIssueKind::unbonded_residue_atoms:
        return "unbonded_residue_atoms";
    }
    throw std::invalid_argument{"unknown import issue kind"};
}

auto aromatic_bonds_issue(std::vector<std::size_t> bonds) -> ImportIssue {
    const auto count = bonds.size();
    const auto first = bonds.front();
    return ImportIssue{
        .kind = ImportIssueKind::aromatic_bonds_converted,
        .message = std::to_string(count) +
                   (count == 1 ? " aromatic bond was" : " aromatic bonds were") +
                   " read as single bonds; methods that use bond orders then give other charges "
                   "than for the Kekulé structure. Write the file with Kekulé bonds to avoid "
                   "this.",
        .line = std::nullopt,
        .location = {.bond_index = first, .bond_indices = std::move(bonds)}};
}

auto make_record(std::vector<core::Atom> atoms, std::vector<core::Bond> bonds,
                 std::vector<core::Conformer> conformers, MoleculeRecordIdentity identity,
                 std::string name, std::vector<ImportIssue> issues,
                 std::optional<MoleculeImportMetadata> import_metadata) -> ImportedMoleculeRecord {
    if (name.empty()) {
        name = identity.record_id.display_string();
    }

    return ImportedMoleculeRecord{.molecule =
                                      core::Molecule{std::move(atoms), std::move(bonds),
                                                     std::move(conformers), std::move(name)},
                                  .identity = std::move(identity),
                                  .caller_atom_ids = std::nullopt,
                                  .import_warnings = std::move(issues),
                                  .import_metadata = std::move(import_metadata)};
}

} // namespace chargefw::adapters
