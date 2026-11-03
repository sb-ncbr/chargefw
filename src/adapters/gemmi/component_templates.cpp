#include "component_templates.h"

#include <algorithm>
#include <array>

namespace chargefw::adapters::gemmi::component_templates {
namespace {

#include "component_templates.inc"

// Alternative Lewis structures follow the wwPDB amino acid protonation variants
// (aa-variants-v1.cif). The CCD templates have the double bond at ND1=CE1 in histidine and at
// CZ=NH2 in arginine.
constexpr std::array<BondTemplate, 2> his_ce1_ne2_double_bonds{{
    {"ND1", "CE1", core::BondOrder::SINGLE},
    {"CE1", "NE2", core::BondOrder::DOUBLE},
}};
constexpr std::array<std::span<const BondTemplate>, 1> his_alternatives{his_ce1_ne2_double_bonds};
constexpr std::array<std::string_view, 5> his_motif_atoms{"CG", "ND1", "CD2", "CE1", "NE2"};
constexpr ProtonationMotif his_motif{.atoms = his_motif_atoms, .alternatives = his_alternatives};

constexpr std::array<BondTemplate, 2> arg_cz_nh1_double_bonds{{
    {"CZ", "NH1", core::BondOrder::DOUBLE},
    {"CZ", "NH2", core::BondOrder::SINGLE},
}};
// Arginine without HE is not a wwPDB variant; force fields name it, for example, AR0.
constexpr std::array<BondTemplate, 2> arg_ne_cz_double_bonds{{
    {"NE", "CZ", core::BondOrder::DOUBLE},
    {"CZ", "NH2", core::BondOrder::SINGLE},
}};
constexpr std::array<std::span<const BondTemplate>, 2> arg_alternatives{arg_cz_nh1_double_bonds,
                                                                        arg_ne_cz_double_bonds};
constexpr std::array<std::string_view, 5> arg_motif_atoms{"CD", "NE", "CZ", "NH1", "NH2"};
constexpr ProtonationMotif arg_motif{.atoms = arg_motif_atoms, .alternatives = arg_alternatives};

struct ComponentTemplateEntry {
    std::string_view component;
    ComponentTemplate component_template;
};

constexpr std::array<ComponentTemplateEntry, 29> templates{{
    {"A", {ComponentKind::nucleotide, a_bonds}},
    {"ALA", {ComponentKind::amino_acid, ala_bonds}},
    {"ARG", {ComponentKind::amino_acid, arg_bonds, &arg_motif}},
    {"ASN", {ComponentKind::amino_acid, asn_bonds}},
    {"ASP", {ComponentKind::amino_acid, asp_bonds}},
    {"C", {ComponentKind::nucleotide, c_bonds}},
    {"CYS", {ComponentKind::amino_acid, cys_bonds}},
    {"DA", {ComponentKind::nucleotide, da_bonds}},
    {"DC", {ComponentKind::nucleotide, dc_bonds}},
    {"DG", {ComponentKind::nucleotide, dg_bonds}},
    {"DT", {ComponentKind::nucleotide, dt_bonds}},
    {"G", {ComponentKind::nucleotide, g_bonds}},
    {"GLN", {ComponentKind::amino_acid, gln_bonds}},
    {"GLU", {ComponentKind::amino_acid, glu_bonds}},
    {"GLY", {ComponentKind::amino_acid, gly_bonds}},
    {"HIS", {ComponentKind::amino_acid, his_bonds, &his_motif}},
    {"HOH", {ComponentKind::water, hoh_bonds}},
    {"ILE", {ComponentKind::amino_acid, ile_bonds}},
    {"LEU", {ComponentKind::amino_acid, leu_bonds}},
    {"LYS", {ComponentKind::amino_acid, lys_bonds}},
    {"MET", {ComponentKind::amino_acid, met_bonds}},
    {"PHE", {ComponentKind::amino_acid, phe_bonds}},
    {"PRO", {ComponentKind::amino_acid, pro_bonds}},
    {"SER", {ComponentKind::amino_acid, ser_bonds}},
    {"THR", {ComponentKind::amino_acid, thr_bonds}},
    {"TRP", {ComponentKind::amino_acid, trp_bonds}},
    {"TYR", {ComponentKind::amino_acid, tyr_bonds}},
    {"U", {ComponentKind::nucleotide, u_bonds}},
    {"VAL", {ComponentKind::amino_acid, val_bonds}},
}};

} // namespace

auto find(const std::string_view component) -> const ComponentTemplate* {
    const auto entry =
        std::ranges::lower_bound(templates, component, {}, &ComponentTemplateEntry::component);
    if (entry != templates.end() && entry->component == component) {
        return std::addressof(entry->component_template);
    }
    return nullptr;
}

} // namespace chargefw::adapters::gemmi::component_templates
