#include <chargefw/adapters/molecule_record.h>

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
    case ImportIssueKind::v2000_property_ignored:
        return "v2000_property_ignored";
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

} // namespace chargefw::adapters
