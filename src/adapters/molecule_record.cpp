#include <chargefw/adapters/molecule_record.h>

#include <stdexcept>
#include <string_view>

namespace chargefw::adapters {

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
