#include <chargefw/io/mmcif_output.h>

#include "io/component.h"
#include "io/gemmi/component_templates.h"
#include "io/gemmi/mmcif_attachment.h"
#include "io/gemmi/mmcif_labels.h"
#include "io/gemmi/pdb_conversion.h"
#include "io/output_helpers.h"
#include "io/record_assignments.h"

#include <chargefw/core/issue_location.h>
#include <chargefw/core/periodic_table.h>
#include <chargefw/io/output_error.h>
#include <chargefw/version.h>

#include <gemmi/cif.hpp>
#include <gemmi/elem.hpp>
#include <gemmi/numb.hpp>
#include <gemmi/to_cif.hpp>
#include <gemmi/util.hpp>

#include <algorithm>
#include <cctype>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <format>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <ostream>
#include <ranges>
#include <set>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace chargefw::io {
namespace {

constexpr auto metadata_category = "_sb_ncbr_partial_atomic_charges_meta.";
constexpr auto charges_category = "_sb_ncbr_partial_atomic_charges.";
constexpr auto dictionary_name = "mmcif_charges_v11.dic";
constexpr auto dictionary_version = "1.1";
constexpr auto dictionary_location =
    "https://sb-ncbr.github.io/charges-schema/schemas/mmcif_charges_v11.dic";
[[nodiscard]] auto atom_site_columns() -> const std::vector<std::string>& {
    static const auto columns = std::vector<std::string>{"id",
                                                         "type_symbol",
                                                         "label_atom_id",
                                                         "label_alt_id",
                                                         "label_comp_id",
                                                         "label_asym_id",
                                                         "label_entity_id",
                                                         "label_seq_id",
                                                         "pdbx_PDB_ins_code",
                                                         "Cartn_x",
                                                         "Cartn_y",
                                                         "Cartn_z",
                                                         "occupancy",
                                                         "B_iso_or_equiv",
                                                         "pdbx_formal_charge",
                                                         "auth_seq_id",
                                                         "auth_comp_id",
                                                         "auth_asym_id",
                                                         "auth_atom_id",
                                                         "pdbx_PDB_model_num"};
    return columns;
}

[[nodiscard]] auto metadata_columns() -> const std::vector<std::string>& {
    static const auto columns = std::vector<std::string>{
        "id", "type", "method", "parameter_set", "software_name", "software_version"};
    return columns;
}

[[nodiscard]] auto charge_columns() -> const std::vector<std::string>& {
    static const auto columns = std::vector<std::string>{"type_id", "atom_id", "charge"};
    return columns;
}

struct BlockMapping {
    std::vector<std::vector<std::string>> atom_site_ids;
};

auto ensure_dictionary(::gemmi::cif::Block& block) -> void;

[[nodiscard]] auto quote(const std::string_view value) -> std::string {
    return ::gemmi::cif::quote(std::string{value});
}

[[nodiscard]] auto block_name(const ImportedMoleculeRecord& record, const std::size_t index)
    -> std::string {
    auto candidate = detail::record_display_name(record, index);
    for (auto& character : candidate) {
        if (std::isalnum(static_cast<unsigned char>(character)) == 0 && character != '_' &&
            character != '-' && character != '.') {
            character = '_';
        }
    }
    return candidate;
}

// Block names must be unique, so repeated names get _2, _3, and so on. CIF block names are
// case-insensitive, so names are compared in lowercase and keep their spelling. Remembering used
// names and the next suffix per name keeps this linear in the number of records, which often share
// a name.
class BlockNames {
  public:
    [[nodiscard]] auto unique(std::string candidate) -> std::string {
        auto key = ::gemmi::to_lower(candidate);
        if (used_.insert(key).second) {
            return candidate;
        }
        auto& suffix = next_suffix_.try_emplace(std::move(key), 2).first->second;
        while (true) {
            auto name = candidate + "_" + std::to_string(suffix++);
            if (used_.insert(::gemmi::to_lower(name)).second) {
                return name;
            }
        }
    }

  private:
    std::unordered_set<std::string> used_;
    std::unordered_map<std::string, std::size_t> next_suffix_;
};

[[nodiscard]] auto atom_ids(const core::Molecule& molecule) -> std::vector<std::string> {
    std::vector<std::string> result;
    result.reserve(molecule.atom_count());
    std::unordered_set<std::string> used;
    for (std::size_t index = 0; index < molecule.atom_count(); ++index) {
        auto candidate = std::string{molecule.atom(index).name()};
        const auto valid =
            !candidate.empty() &&
            std::ranges::all_of(candidate, [](const unsigned char character) -> bool {
                return std::isalnum(character) != 0 || character == '_' || character == '-';
            });
        if (!valid || used.contains(candidate)) {
            candidate = std::string{core::element_symbol(molecule.atom(index).atomic_number())} +
                        std::to_string(index + 1);
        }
        while (used.contains(candidate)) {
            candidate += "_" + std::to_string(index + 1);
        }
        used.insert(candidate);
        result.push_back(std::move(candidate));
    }
    return result;
}

struct OutputAssignments {
    const charges::ChargeAssignment* molecule = nullptr;
    std::vector<const charges::ChargeAssignment*> conformers;
};

[[nodiscard]] auto assignments_for_result(const ChargeCalculationResult& result,
                                          const std::size_t molecule_index) -> OutputAssignments {
    const auto& molecule = result.inputs()[molecule_index].molecule;
    auto assignments = OutputAssignments{
        .molecule = nullptr,
        .conformers = std::vector<const charges::ChargeAssignment*>(molecule.conformer_count())};
    for (const auto& assignment : record_assignments(result, molecule_index)) {
        if (assignment.target.conformer_index.has_value()) {
            assignments.conformers[*assignment.target.conformer_index] = std::addressof(assignment);
        } else {
            assignments.molecule = std::addressof(assignment);
        }
    }
    return assignments;
}

[[nodiscard]] auto canonical_charge_site_id(const std::string_view value) -> bool {
    int id = 0;
    const auto [end, error] = std::from_chars(value.data(), value.data() + value.size(), id);
    return error == std::errc{} && end == value.data() + value.size() && id > 0 &&
           std::to_string(id) == value;
}

// Where a refusal of one atom of a conformer is; the molecule is added by rethrow_for_molecule().
[[nodiscard]] auto atom_location(const std::size_t conformer_index, const std::size_t atom_index)
    -> core::IssueLocation {
    return core::IssueLocation{
        .conformer_index = conformer_index, .atom_index = atom_index, .atom_indices = {atom_index}};
}

// Where a refusal of two atoms is, in index order.
[[nodiscard]] auto atom_pair(const std::size_t first, const std::size_t second)
    -> core::IssueLocation {
    const auto [lower, upper] = std::minmax(first, second);
    return core::IssueLocation{
        .atom_index = lower, .other_atom_index = upper, .atom_indices = {lower, upper}};
}

// Called from a catch block around the work on one input molecule: a refusal that does not name
// its molecule gets it.
[[noreturn]] auto rethrow_for_molecule(const std::size_t molecule_index) -> void {
    try {
        throw;
    } catch (const OutputError& error) {
        if (error.molecule_index().has_value()) {
            throw;
        }
        throw OutputError{error.kind(), error.reason(), molecule_index, error.location()};
    }
}

[[nodiscard]] auto attached_mapping(::gemmi::cif::Block& block,
                                    const ImportedMoleculeRecord& record) -> BlockMapping {
    if (!record.import_metadata.has_value() ||
        record.import_metadata->format != MolecularSourceFormat::mmcif) {
        throw OutputError{
            OutputErrorKind::unsupported_source,
            "charges can be added only to the mmCIF document the molecules were read from"};
    }
    if (record.identity.record_id.empty() ||
        record.identity.record_id.display_string() != block.name) {
        throw OutputError{OutputErrorKind::target_mismatch,
                          "Gemmi target block identity does not match the calculation input"};
    }
    auto atom_sites = gemmi::mmcif_labels::source_atom_sites(block);
    auto type_symbols = block.find("_atom_site.", {"type_symbol"});
    auto identity_sites =
        block.find("_atom_site.", {"Cartn_x", "Cartn_y", "Cartn_z", "?pdbx_formal_charge"});
    if (atom_sites.length() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw OutputError{OutputErrorKind::unrepresentable,
                          "Gemmi target contains too many atom sites"};
    }
    if (type_symbols.length() != atom_sites.length() ||
        identity_sites.length() != atom_sites.length()) {
        throw OutputError{OutputErrorKind::target_mismatch,
                          "Gemmi target site mapping does not match the calculation input"};
    }
    const auto& metadata = *record.import_metadata;
    auto result = BlockMapping{};
    result.atom_site_ids.reserve(metadata.conformers.size());
    auto mapped_ids = std::unordered_set<std::string>{};
    for (std::size_t conformer_index = 0; conformer_index < metadata.conformers.size();
         ++conformer_index) {
        const auto& conformer = metadata.conformers[conformer_index];
        const auto model_id = conformer.id.value_or("1");
        auto& ids = result.atom_site_ids.emplace_back();
        ids.reserve(conformer.sites.size());
        for (std::size_t atom_index = 0; atom_index < conformer.sites.size(); ++atom_index) {
            const auto& site = conformer.sites[atom_index];
            if (!site.id.has_value() || site.position >= atom_sites.length()) {
                throw OutputError{OutputErrorKind::target_mismatch,
                                  "Gemmi target site mapping does not match the calculation input",
                                  std::nullopt, atom_location(conformer_index, atom_index)};
            }
            const auto row_index = static_cast<int>(site.position);
            auto target_id = ::gemmi::cif::as_string(atom_sites[row_index][0]);
            if (target_id != *site.id || !mapped_ids.insert(target_id).second) {
                throw OutputError{OutputErrorKind::target_mismatch,
                                  "Gemmi target site mapping does not match the calculation input",
                                  std::nullopt, atom_location(conformer_index, atom_index)};
            }
            if (!canonical_charge_site_id(target_id)) {
                throw OutputError{
                    OutputErrorKind::unrepresentable,
                    "Gemmi target atom IDs are not representable by the charge dictionary",
                    std::nullopt, atom_location(conformer_index, atom_index)};
            }
            const auto& target_site = atom_sites[row_index];
            // Compared as elements, as the reader reads them: wwPDB files write SE for Se.
            if (::gemmi::Element(::gemmi::cif::as_string(type_symbols[row_index][0]))
                    .atomic_number() != record.molecule.atom(atom_index).atomic_number()) {
                throw OutputError{OutputErrorKind::target_mismatch,
                                  "Gemmi target site mapping does not match the calculation input",
                                  std::nullopt, atom_location(conformer_index, atom_index)};
            }
            if (site.structural_labels == nullptr ||
                gemmi::mmcif_labels::decode_source_labels(target_site) != *site.structural_labels) {
                throw OutputError{OutputErrorKind::target_mismatch,
                                  "Gemmi target site mapping does not match the calculation input",
                                  std::nullopt, atom_location(conformer_index, atom_index)};
            }
            if (target_site.has(gemmi::mmcif_labels::model_id_column) &&
                ::gemmi::cif::as_string(target_site[gemmi::mmcif_labels::model_id_column]) !=
                    model_id) {
                throw OutputError{OutputErrorKind::target_mismatch,
                                  "Gemmi target model mapping does not match the calculation input",
                                  std::nullopt, atom_location(conformer_index, atom_index)};
            }
            const auto& position = record.molecule.conformer(conformer_index)[atom_index];
            const auto& identity_site = identity_sites[row_index];
            // As the reader reads it, and zero when the column is absent. Text that is not a
            // formal charge is a changed target.
            const auto formal_charge = identity_site.has(3)
                                           ? gemmi::mmcif_labels::formal_charge(identity_site[3])
                                           : std::optional<int>{0};
            if (::gemmi::cif::as_number(identity_site[0]) != position.x ||
                ::gemmi::cif::as_number(identity_site[1]) != position.y ||
                ::gemmi::cif::as_number(identity_site[2]) != position.z ||
                formal_charge != record.molecule.atom(atom_index).formal_charge()) {
                throw OutputError{OutputErrorKind::target_mismatch,
                                  "Gemmi target site mapping does not match the calculation input",
                                  std::nullopt, atom_location(conformer_index, atom_index)};
            }
            ids.push_back(std::move(target_id));
        }
    }
    return result;
}

[[nodiscard]] auto structural_labels(const ImportedMoleculeRecord& record,
                                     const std::size_t conformer_index,
                                     const std::size_t atom_index)
    -> const SourceStructuralLabels* {
    if (!record.import_metadata.has_value()) {
        return nullptr;
    }
    const auto& reference = record.import_metadata->conformers[conformer_index].sites[atom_index];
    return reference.structural_labels.get();
}

struct OutputSiteLabels {
    std::string label_atom;
    std::string author_atom;
    std::string label_residue;
    std::string author_residue;
    std::string label_chain;
    std::string author_chain;
    std::string label_sequence;
    std::string author_sequence;
    std::string entity;
    std::string alternate_location;
    std::string insertion_code;
};

[[nodiscard]] auto output_site_labels(const ImportedMoleculeRecord& record,
                                      const std::vector<std::string>& generated_atom_ids,
                                      const std::size_t conformer_index,
                                      const std::size_t atom_index) -> OutputSiteLabels {
    const auto* labels = structural_labels(record, conformer_index, atom_index);
    const auto& fallback_atom = generated_atom_ids[atom_index];
    if (labels == nullptr) {
        return {.label_atom = fallback_atom,
                .author_atom = fallback_atom,
                .label_residue = "UNL",
                .author_residue = "UNL",
                .label_chain = "A",
                .author_chain = "A",
                .label_sequence = ".",
                .author_sequence = "1",
                .entity = "1",
                .alternate_location = ".",
                .insertion_code = "."};
    }
    const auto known_or = [](const std::optional<std::string>& value,
                             const std::string_view fallback) {
        return std::string{detail::known_label(value).value_or(fallback)};
    };
    const auto label_or = [](const std::optional<std::string>& preferred,
                             const std::optional<std::string>& alternate,
                             const std::string_view fallback) {
        return std::string{detail::first_known_label(preferred, alternate, fallback)};
    };
    return {.label_atom = label_or(labels->label.atom, labels->author.atom, fallback_atom),
            .author_atom = label_or(labels->author.atom, labels->label.atom, fallback_atom),
            .label_residue = label_or(labels->label.residue, labels->author.residue, "UNL"),
            .author_residue = label_or(labels->author.residue, labels->label.residue, "UNL"),
            .label_chain = label_or(labels->label.chain, labels->author.chain, "A"),
            .author_chain = label_or(labels->author.chain, labels->label.chain, "A"),
            .label_sequence = known_or(labels->label.sequence, "."),
            .author_sequence = known_or(labels->author.sequence, "."),
            .entity = known_or(labels->entity, "1"),
            .alternate_location = known_or(labels->alternate_location, "."),
            .insertion_code = known_or(labels->insertion_code, ".")};
}

// The bond order codes of _chem_comp_bond.value_order and _struct_conn.pdbx_value_order.
[[nodiscard]] auto mmcif_bond_order(const core::BondOrder order) -> std::string_view {
    switch (order) {
    case core::BondOrder::SINGLE:
        return "sing";
    case core::BondOrder::DOUBLE:
        return "doub";
    case core::BondOrder::TRIPLE:
        return "trip";
    }
    throw std::logic_error{"unknown bond order"};
}

// The residue instance of each atom. Structural imports record the selected residues; any other
// record is one residue.
[[nodiscard]] auto residue_instances(const ImportedMoleculeRecord& record)
    -> std::vector<std::size_t> {
    const auto atom_count = record.molecule.atom_count();
    if (!record.import_metadata.has_value() || record.import_metadata->components.empty()) {
        return std::vector<std::size_t>(atom_count, 0);
    }
    constexpr auto missing = std::numeric_limits<std::size_t>::max();
    auto result = std::vector<std::size_t>(atom_count, missing);
    const auto& components = record.import_metadata->components;
    for (std::size_t component_index = 0; component_index < components.size(); ++component_index) {
        for (const auto atom_index : components[component_index].atom_indices) {
            if (atom_index >= atom_count || result[atom_index] != missing) {
                throw std::invalid_argument{
                    "mmCIF bond output requires each atom in exactly one residue"};
            }
            result[atom_index] = component_index;
        }
    }
    if (std::ranges::contains(result, missing)) {
        throw std::invalid_argument{"mmCIF bond output requires each atom in exactly one residue"};
    }
    return result;
}

[[nodiscard]] auto residue_description(const OutputSiteLabels& labels) -> std::string {
    return "residue " + detail::residue_description(labels.label_residue, labels.author_chain,
                                                    labels.author_sequence, labels.insertion_code);
}

// Two atom names in sorted order, so that a bond has one key whichever atom comes first.
using AtomNamePair = std::pair<std::string, std::string>;

[[nodiscard]] auto bond_key(std::string first, std::string second) -> AtomNamePair {
    if (second < first) {
        std::swap(first, second);
    }
    return {std::move(first), std::move(second)};
}

struct ResidueBonds {
    std::string component;
    std::string description;
    std::map<std::string, std::size_t, std::less<>> atoms;
    std::map<AtomNamePair, core::BondOrder> bonds;
};

struct ComponentBond {
    std::string component;
    std::string first_atom;
    std::string second_atom;
    core::BondOrder order;
    bool delocalized = false;
};

struct RecordBonds {
    std::vector<ComponentBond> component_bonds;
    std::vector<const core::Bond*> connections;
};

// Whether every residue with the bond's name has the order that ChargeFW's template gives the bond
// for the residue's hydrogens and formal charges, as for the double bond of the histidine ring.
[[nodiscard]] auto follows_protonation(const core::Molecule& molecule,
                                       const std::vector<ResidueBonds>& residues,
                                       const ComponentBond& bond) -> bool {
    const auto* const component_template = gemmi::component_templates::find(bond.component);
    if (component_template == nullptr ||
        !gemmi::component_templates::varies_with_protonation(*component_template, bond.first_atom,
                                                             bond.second_atom)) {
        return false;
    }
    const auto key = bond_key(bond.first_atom, bond.second_atom);
    return std::ranges::all_of(residues, [&](const ResidueBonds& residue) {
        const auto order = residue.bonds.find(key);
        if (residue.component != bond.component || order == residue.bonds.end()) {
            return true;
        }
        const auto atom = [&](const std::string_view name)
            -> std::optional<gemmi::component_templates::MotifAtom> {
            const auto found = residue.atoms.find(name);
            if (found == residue.atoms.end()) {
                return std::nullopt;
            }
            const auto& source = molecule.atom(found->second);
            return gemmi::component_templates::MotifAtom{.atomic_number = source.atomic_number(),
                                                         .formal_charge = source.formal_charge()};
        };
        const auto protonation =
            gemmi::component_templates::protonation_bonds(*component_template, atom);
        return gemmi::component_templates::bond_order(*component_template, protonation,
                                                      bond.first_atom,
                                                      bond.second_atom) == order->second;
    });
}

// Splits the bonds into component bonds, written per residue name, and connections between
// residues. A component bond applies to every residue with that name, so all such residues must
// have the same bonds between the atoms they share. Residues may differ only in the order of a bond
// that follows their protonation state; that bond is written as delocalized, and ChargeFW's
// templates give each residue its order again when the file is read.
[[nodiscard]] auto record_bonds(const core::Molecule& molecule,
                                const std::vector<OutputSiteLabels>& labels,
                                const std::vector<std::size_t>& residue_of) -> RecordBonds {
    auto result = RecordBonds{};
    auto residues = std::vector<ResidueBonds>(std::ranges::max(residue_of) + 1);
    for (std::size_t atom_index = 0; atom_index < molecule.atom_count(); ++atom_index) {
        const auto& site = labels[atom_index];
        auto& residue = residues[residue_of[atom_index]];
        if (residue.atoms.empty()) {
            residue.component = site.label_residue;
            residue.description = residue_description(site);
        }
        if (const auto [named, added] = residue.atoms.emplace(site.label_atom, atom_index);
            !added) {
            throw OutputError{
                OutputErrorKind::unrepresentable,
                std::format("mmCIF bond output requires unique atom names, but {} has atom '{}' "
                            "more than once",
                            residue.description, site.label_atom),
                std::nullopt, atom_pair(named->second, atom_index)};
        }
    }

    // The index of each component bond, and the first residue whose order differs from it.
    auto definitions = std::map<std::string, std::map<AtomNamePair, std::size_t>>{};
    auto differing = std::map<std::size_t, std::string>{};
    for (const auto& bond : molecule.bonds()) {
        const auto residue_index = residue_of[bond.first_atom_index()];
        if (residue_index != residue_of[bond.second_atom_index()]) {
            result.connections.push_back(std::addressof(bond));
            continue;
        }
        auto& residue = residues[residue_index];
        const auto& first = labels[bond.first_atom_index()].label_atom;
        const auto& second = labels[bond.second_atom_index()].label_atom;
        const auto key = bond_key(first, second);
        residue.bonds.emplace(key, bond.order());
        const auto [definition, inserted] =
            definitions[residue.component].try_emplace(key, result.component_bonds.size());
        if (inserted) {
            result.component_bonds.push_back({.component = residue.component,
                                              .first_atom = first,
                                              .second_atom = second,
                                              .order = bond.order()});
        } else if (result.component_bonds[definition->second].order != bond.order()) {
            differing.try_emplace(definition->second, residue.description);
        }
    }

    for (const auto& [index, description] : differing) {
        auto& bond = result.component_bonds[index];
        if (!follows_protonation(molecule, residues, bond)) {
            throw OutputError{
                OutputErrorKind::unrepresentable,
                std::format("mmCIF bond output cannot write bond {}-{} of {}: another {} residue "
                            "has a different bond order",
                            bond.first_atom, bond.second_atom, description, bond.component)};
        }
        bond.delocalized = true;
    }

    // Rows of a residue name follow its template where it has one, because reading the file back
    // puts the bonds the file states before those the template adds, such as delocalized ones.
    auto component_rank = std::map<std::string, std::size_t, std::less<>>{};
    for (const auto& bond : result.component_bonds) {
        component_rank.try_emplace(bond.component, component_rank.size());
    }
    const auto template_position = [](const ComponentBond& bond) {
        constexpr auto untemplated = std::numeric_limits<std::size_t>::max();
        const auto* const component_template = gemmi::component_templates::find(bond.component);
        if (component_template == nullptr) {
            return untemplated;
        }
        const auto found =
            std::ranges::find_if(component_template->bonds, [&bond](const auto& row) {
                return bond_key(std::string{row.first}, std::string{row.second}) ==
                       bond_key(bond.first_atom, bond.second_atom);
            });
        return found == component_template->bonds.end()
                   ? untemplated
                   : static_cast<std::size_t>(found - component_template->bonds.begin());
    };
    std::ranges::stable_sort(result.component_bonds, {}, [&](const ComponentBond& bond) {
        return std::pair{component_rank.find(bond.component)->second, template_position(bond)};
    });

    // A residue without a bond that its component defines would gain that bond in readers.
    for (const auto& residue : residues) {
        for (const auto& key : std::views::keys(definitions[residue.component])) {
            if (residue.atoms.contains(key.first) && residue.atoms.contains(key.second) &&
                !residue.bonds.contains(key)) {
                throw OutputError{
                    OutputErrorKind::unrepresentable,
                    std::format("mmCIF bond output cannot write bond {}-{} of {}: {} has both "
                                "atoms but no bond between them",
                                key.first, key.second, residue.component, residue.description),
                    std::nullopt,
                    atom_pair(residue.atoms.find(key.first)->second,
                              residue.atoms.find(key.second)->second)};
            }
        }
    }
    return result;
}

// No _chem_comp rows are written: _chem_comp requires each component's type, which ChargeFW does
// not know. Omitting the category lets readers that support residue-name fallback use it.
auto append_component_bonds(::gemmi::cif::Block& block, const RecordBonds& bonds) -> void {
    if (bonds.component_bonds.empty()) {
        return;
    }
    const auto columns = std::vector<std::string>{"comp_id", "atom_id_1", "atom_id_2",
                                                  "value_order", "pdbx_ordinal"};
    block.init_loop("_chem_comp_bond.", columns);
    auto rows = block.find("_chem_comp_bond.", columns);
    auto ordinal = std::size_t{1};
    for (const auto& bond : bonds.component_bonds) {
        rows.append_row({quote(bond.component), quote(bond.first_atom), quote(bond.second_atom),
                         std::string{bond.delocalized ? "delo" : mmcif_bond_order(bond.order)},
                         std::to_string(ordinal++)});
    }
}

// Appends the _struct_conn partner columns of one atom, in the order of connection_columns().
// A refusal is located at the two atoms of the bond.
auto append_connection_partner(std::vector<std::string>& row, const OutputSiteLabels& site,
                               const core::IssueLocation& bond_atoms) -> void {
    // Readers find a partner by author chain and number, or by label chain and number.
    if (site.author_sequence == "." && site.label_sequence == ".") {
        throw OutputError{OutputErrorKind::unrepresentable,
                          "mmCIF bond output cannot write a bond to " + residue_description(site) +
                              ": it has no residue number",
                          std::nullopt, bond_atoms};
    }
    // Readers such as Gemmi name atoms by their author names but look up partners by the label
    // names, so a partner whose names differ would not be found.
    if (site.label_atom != site.author_atom || site.label_residue != site.author_residue) {
        throw OutputError{
            OutputErrorKind::unrepresentable,
            std::format("mmCIF bond output cannot write a bond to atom {} of {}: its author and "
                        "label atom or residue names differ",
                        site.label_atom, residue_description(site)),
            std::nullopt, bond_atoms};
    }
    const auto optional = [](const std::string& value) {
        return value == "." ? std::string{"."} : quote(value);
    };
    row.insert(row.end(),
               {quote(site.label_chain), quote(site.label_residue), optional(site.label_sequence),
                quote(site.label_atom), optional(site.insertion_code), quote(site.author_chain),
                optional(site.author_sequence), "1_555"});
}

[[nodiscard]] auto connection_columns() -> std::vector<std::string> {
    auto columns = std::vector<std::string>{"id", "conn_type_id"};
    for (const auto* partner : {"1", "2"}) {
        const auto ptnr = std::string{"ptnr"} + partner;
        const auto pdbx_ptnr = std::string{"pdbx_ptnr"} + partner;
        columns.insert(columns.end(),
                       {ptnr + "_label_asym_id", ptnr + "_label_comp_id", ptnr + "_label_seq_id",
                        ptnr + "_label_atom_id", pdbx_ptnr + "_PDB_ins_code",
                        ptnr + "_auth_asym_id", ptnr + "_auth_seq_id", ptnr + "_symmetry"});
    }
    columns.emplace_back("pdbx_value_order");
    return columns;
}

auto append_connections(::gemmi::cif::Block& block, const core::Molecule& molecule,
                        const std::vector<OutputSiteLabels>& labels,
                        const std::vector<const core::Bond*>& connections) -> void {
    if (connections.empty()) {
        return;
    }
    const auto columns = connection_columns();
    block.init_loop("_struct_conn.", columns);
    auto rows = block.find("_struct_conn.", columns);
    auto next_covalent_id = std::size_t{1};
    auto next_disulfide_id = std::size_t{1};
    for (const auto* bond : connections) {
        constexpr auto sulfur = 16;
        const auto disulfide = molecule.atom(bond->first_atom_index()).atomic_number() == sulfur &&
                               molecule.atom(bond->second_atom_index()).atomic_number() == sulfur;
        const auto type = std::string{disulfide ? "disulf" : "covale"};
        auto row = std::vector<std::string>{
            type + std::to_string(disulfide ? next_disulfide_id++ : next_covalent_id++), type};
        const auto bond_atoms = atom_pair(bond->first_atom_index(), bond->second_atom_index());
        append_connection_partner(row, labels[bond->first_atom_index()], bond_atoms);
        append_connection_partner(row, labels[bond->second_atom_index()], bond_atoms);
        row.emplace_back(mmcif_bond_order(bond->order()));
        rows.append_row(row);
    }
}

// mmCIF describes bonds within a residue in _chem_comp_bond, by residue and atom name, and bonds
// between residues in _struct_conn. Bonds are topology and use the first conformer's labels.
auto append_bonds(::gemmi::cif::Block& block, const ImportedMoleculeRecord& record,
                  const std::vector<std::string>& generated_atom_ids) -> void {
    const auto& molecule = record.molecule;
    if (molecule.bond_count() == 0) {
        return;
    }
    auto labels = std::vector<OutputSiteLabels>{};
    labels.reserve(molecule.atom_count());
    for (std::size_t atom_index = 0; atom_index < molecule.atom_count(); ++atom_index) {
        labels.push_back(output_site_labels(record, generated_atom_ids, 0, atom_index));
    }
    const auto bonds = record_bonds(molecule, labels, residue_instances(record));
    append_component_bonds(block, bonds);
    append_connections(block, molecule, labels, bonds.connections);
}

// What the format can hold of a result that require_results() accepted.
auto validate_result_output(const ChargeCalculationResult& result) -> void {
    for (std::size_t molecule_index = 0; molecule_index < result.inputs().size();
         ++molecule_index) {
        const auto& molecule = result.inputs()[molecule_index].molecule;
        if (molecule.atom_count() == 0) {
            throw OutputError{OutputErrorKind::unsupported_shape,
                              "mmCIF output requires at least one atom", molecule_index};
        }
        if (molecule.conformer_count() == 0) {
            throw OutputError{OutputErrorKind::missing_coordinates,
                              "mmCIF output requires coordinates", molecule_index};
        }
        if (molecule.conformer_count() >
            static_cast<std::size_t>(std::numeric_limits<int>::max()) / molecule.atom_count()) {
            throw OutputError{OutputErrorKind::unrepresentable,
                              "mmCIF output contains too many atom sites", molecule_index};
        }
        for (std::size_t conformer_index = 0; conformer_index < molecule.conformer_count();
             ++conformer_index) {
            const auto& conformer = molecule.conformer(conformer_index);
            for (std::size_t atom_index = 0; atom_index < molecule.atom_count(); ++atom_index) {
                const auto& position = conformer[atom_index];
                if (!std::isfinite(position.x) || !std::isfinite(position.y) ||
                    !std::isfinite(position.z)) {
                    throw OutputError{OutputErrorKind::unrepresentable,
                                      "mmCIF output coordinates must be finite", molecule_index,
                                      atom_location(conformer_index, atom_index)};
                }
            }
        }
    }
    for (const auto& assignment : result.charges()->assignments()) {
        const auto values = assignment.charges.values();
        for (std::size_t atom_index = 0; atom_index < values.size(); ++atom_index) {
            if (values[atom_index] < -10.0 || values[atom_index] > 10.0) {
                throw OutputError{
                    OutputErrorKind::unrepresentable,
                    "mmCIF charge is outside the dictionary range [-10, 10]",
                    assignment.target.molecule_index,
                    core::IssueLocation{.conformer_index = assignment.target.conformer_index,
                                        .atom_index = atom_index,
                                        .atom_indices = {atom_index}}};
            }
        }
    }
}

using ResultSpan = std::span<const ChargeCalculationResult* const>;

// Several charge sets share one file only when they were calculated for the same molecules and
// can be told apart; the file records method and parameter set but not options or execution.
// Results that can be written at all: successful, and with molecules.
auto require_results(const ResultSpan results) -> void {
    if (results.empty()) {
        throw OutputError{OutputErrorKind::unsupported_shape,
                          "mmCIF output requires at least one calculation result"};
    }
    for (const auto* result : results) {
        if (result == nullptr) {
            throw std::invalid_argument{"mmCIF output requires calculation results, not null"};
        }
        if (result->charges() == nullptr) {
            throw OutputError{OutputErrorKind::no_charges,
                              "mmCIF output requires a successful calculation"};
        }
        if (result->inputs().empty()) {
            throw OutputError{OutputErrorKind::unsupported_shape,
                              "mmCIF output requires at least one molecule"};
        }
    }
}

auto validate_results_output(const ResultSpan results) -> void {
    require_results(results);
    auto charge_sets = std::set<std::pair<std::string, std::optional<std::string>>>{};
    for (const auto* result : results) {
        validate_result_output(*result);
        if (result != results.front() &&
            !std::ranges::equal(result->inputs(), results.front()->inputs())) {
            throw std::invalid_argument{
                "mmCIF output requires charge sets calculated for the same molecules"};
        }
        const auto& charge_set = *result->charges();
        const auto parameter_set = charge_set.parameter_set_id().transform(
            [](const std::string_view value) { return std::string{value}; });
        if (!charge_sets.emplace(std::string{charge_set.method_id()}, parameter_set).second) {
            throw OutputError{
                OutputErrorKind::unrepresentable,
                std::format("mmCIF output cannot hold two charge sets of method {} with parameter "
                            "set {}: the file cannot distinguish their options or execution",
                            charge_set.method_id(), parameter_set.value_or("none"))};
        }
    }
}

auto append_charge_metadata(::gemmi::cif::Table& metadata, const std::string_view id,
                            const charges::ChargeSet& charge_set) -> void {
    const auto parameter_set = charge_set.parameter_set_id();
    metadata.append_row({std::string{id}, "empirical", quote(charge_set.method_id()),
                         parameter_set.has_value() ? quote(*parameter_set) : ".", quote("ChargeFW"),
                         quote(CHARGEFW_VERSION_STRING)});
}

// Checks the charge sets a block already has, for adding the given ones, and returns their IDs.
// Both categories become loops, and the metadata gets the parameter_set column if it lacks it.
[[nodiscard]] auto existing_charge_sets(::gemmi::cif::Block& block, const ResultSpan results)
    -> std::set<std::string> {
    const auto refuse = [](const std::string_view reason) {
        return OutputError{OutputErrorKind::existing_charges,
                           std::format("charge sets cannot be appended: {}", reason)};
    };
    // Rows under another version of the charges dictionary are not reinterpreted.
    for (const auto row : block.find("_audit_conform.", {"dict_name"})) {
        if (const auto name = row.str(0);
            name.starts_with("mmcif_charges_") && name != dictionary_name) {
            throw refuse("the existing charges are declared to follow " + name);
        }
    }
    for (const auto* const category : {metadata_category, charges_category}) {
        if (auto any = block.find_mmcif_category(category); any.ok()) {
            any.ensure_loop();
        }
    }
    auto metadata = block.find(metadata_category, {"id", "type", "method", "?parameter_set",
                                                   "software_name", "software_version"});
    auto charge_rows = block.find(charges_category, charge_columns());
    if (!metadata.ok() || !charge_rows.ok()) {
        throw refuse("the existing charges do not follow dictionary 1.1");
    }
    auto ids = std::set<std::string>{};
    for (const auto row : metadata) {
        if (!ids.insert(row.str(0)).second) {
            throw refuse("the existing charge sets have repeated IDs");
        }
        const auto parameter_set =
            row.has2(3) ? std::optional<std::string>{row.str(3)} : std::nullopt;
        for (const auto* const result : results) {
            const auto& added = *result->charges();
            if (row.str(2) == added.method_id() && parameter_set == added.parameter_set_id()) {
                throw refuse(std::format("there is a charge set of method {} with parameter set {}",
                                         added.method_id(), parameter_set.value_or("none")));
            }
        }
    }
    // Every existing charge belongs to a set and to an atom site of the block, once.
    auto sites = std::unordered_set<std::string>{};
    for (const auto row : block.find("_atom_site.", {"id"})) {
        sites.insert(row.str(0));
    }
    auto charged = std::set<std::pair<std::string, std::string>>{};
    for (const auto row : charge_rows) {
        if (!ids.contains(row.str(0)) || !sites.contains(row.str(1))) {
            throw refuse("an existing charge refers to an unknown charge set or atom site");
        }
        if (!charged.emplace(row.str(0), row.str(1)).second) {
            throw refuse("an atom site has two charges in one existing charge set");
        }
    }
    if (!metadata.has_column(3)) {
        metadata.get_loop()->add_columns({std::string{metadata_category} + "parameter_set"}, ".");
    }
    return ids;
}

// The number that new charge set IDs follow: the largest existing ID that is a number.
[[nodiscard]] auto largest_number(const std::set<std::string>& ids) -> std::uint64_t {
    auto largest = std::uint64_t{0};
    for (const auto& id : ids) {
        auto value = std::uint64_t{0};
        const auto [end, error] = std::from_chars(id.data(), id.data() + id.size(), value);
        if (error == std::errc{} && end == id.data() + id.size()) {
            largest = std::max(largest, value);
        }
    }
    if (largest > std::numeric_limits<std::uint64_t>::max() - 1'000'000) {
        throw OutputError{
            OutputErrorKind::existing_charges,
            "charge sets cannot be appended: the existing charge set IDs are too large"};
    }
    return largest;
}

// Adds one metadata row per charge set and its charges for one record. The sets replace those of
// the block and get IDs from 1, or with keep_existing follow them with the next numbers. Each
// covers every conformer; site_ids[conformer][atom] are the record's _atom_site.id values.
auto append_charge_sets(::gemmi::cif::Block& block, const ResultSpan results,
                        const std::size_t record_index,
                        const std::vector<std::vector<std::string>>& site_ids,
                        const bool keep_existing = false) -> void {
    auto last_id = std::uint64_t{0};
    if (keep_existing && (block.has_mmcif_category(metadata_category) ||
                          block.has_mmcif_category(charges_category))) {
        last_id = largest_number(existing_charge_sets(block, results));
        ensure_dictionary(block);
    } else {
        ensure_dictionary(block);
        block.init_mmcif_loop(metadata_category, metadata_columns());
        block.init_mmcif_loop(charges_category, charge_columns());
    }

    // Category initialization can invalidate tables, so acquire them only after the mutations.
    auto metadata = block.find(metadata_category, metadata_columns());
    auto charge_rows = block.find(charges_category, charge_columns());
    for (std::size_t result_index = 0; result_index < results.size(); ++result_index) {
        const auto& result = *results[result_index];
        const auto id = std::to_string(last_id + result_index + 1);
        append_charge_metadata(metadata, id, *result.charges());
        const auto assignments = assignments_for_result(result, record_index);
        for (std::size_t conformer = 0; conformer < site_ids.size(); ++conformer) {
            const auto* assignment = assignments.molecule != nullptr
                                         ? assignments.molecule
                                         : assignments.conformers.at(conformer);
            if (assignment == nullptr || assignment->charges.size() != site_ids[conformer].size()) {
                throw std::runtime_error{"charge assignment does not match the mmCIF atom sites"};
            }
            for (std::size_t atom = 0; atom < site_ids[conformer].size(); ++atom) {
                charge_rows.append_row({id, quote(site_ids[conformer][atom]),
                                        std::format("{:.4f}", assignment->charges[atom])});
            }
        }
    }
}

auto write_result_block(::gemmi::cif::Block& block, const ResultSpan results,
                        const std::size_t record_index) -> void {
    const auto& record = results.front()->inputs()[record_index];
    const auto& molecule = record.molecule;
    const auto generated_atom_ids = atom_ids(molecule);
    block.set_pair("_entry.id", quote(block.name));
    // Declared here as well as with the charges so that it precedes _atom_site.
    ensure_dictionary(block);

    block.init_loop("_atom_site.", atom_site_columns());
    auto atom_sites = block.find("_atom_site.", atom_site_columns());
    auto site_ids = std::vector<std::vector<std::string>>(molecule.conformer_count());
    auto next_site_id = std::size_t{1};
    for (std::size_t conformer_index = 0; conformer_index < molecule.conformer_count();
         ++conformer_index) {
        const auto model_id = std::to_string(conformer_index + 1);
        for (std::size_t atom_index = 0; atom_index < molecule.atom_count(); ++atom_index) {
            const auto& atom = molecule.atom(atom_index);
            const auto& position = molecule.conformer(conformer_index)[atom_index];
            const auto labels =
                output_site_labels(record, generated_atom_ids, conformer_index, atom_index);
            const auto& site_id =
                site_ids[conformer_index].emplace_back(std::to_string(next_site_id++));
            atom_sites.append_row(
                {site_id,
                 std::string{core::element_symbol(atom.atomic_number())},
                 quote(labels.label_atom),
                 labels.alternate_location == "." ? "." : quote(labels.alternate_location),
                 quote(labels.label_residue),
                 quote(labels.label_chain),
                 quote(labels.entity),
                 labels.label_sequence == "." ? "." : quote(labels.label_sequence),
                 labels.insertion_code == "." ? "." : quote(labels.insertion_code),
                 std::format("{}", position.x),
                 std::format("{}", position.y),
                 std::format("{}", position.z),
                 "?",
                 "?",
                 std::to_string(atom.formal_charge()),
                 labels.author_sequence == "." ? "." : quote(labels.author_sequence),
                 quote(labels.author_residue),
                 quote(labels.author_chain),
                 quote(labels.author_atom),
                 model_id});
        }
    }
    append_charge_sets(block, results, record_index, site_ids);
    append_bonds(block, record, generated_atom_ids);
}

// Declares the charges dictionary once. Declarations of other charges dictionary versions, such as
// mmcif_charges_v10.dic in older files, no longer describe the charges and are removed; other
// dictionaries stay.
auto ensure_dictionary(::gemmi::cif::Block& block) -> void {
    static const auto columns =
        std::vector<std::string>{"dict_name", "dict_version", "dict_location"};
    if (auto category = block.find_mmcif_category("_audit_conform."); category.ok()) {
        // Keep the declarations that are there: add the columns they lack, with unknown values.
        category.ensure_loop();
        auto& loop = *category.get_loop();
        for (const auto& column : columns) {
            if (const auto tag = "_audit_conform." + column; !loop.has_tag(tag)) {
                loop.add_columns({tag}, "?");
            }
        }
    } else {
        block.init_mmcif_loop("_audit_conform.", columns);
    }
    auto table = block.find("_audit_conform.", columns);
    for (auto index = static_cast<int>(table.length()); index-- > 0;) {
        if (::gemmi::cif::as_string(table[index][0]).starts_with("mmcif_charges_")) {
            table.remove_row(index);
        }
    }
    table.append_row({dictionary_name, dictionary_version, dictionary_location});
}

// Adds the charges of validated results to the document and writes it. Record k belongs to the
// block with coordinates number record_index of k.
auto write_charged(std::ostream& output, const ResultSpan results, ::gemmi::cif::Document document,
                   const ExistingCharges existing, const bool every_block) -> void {
    const auto inputs = results.front()->inputs();
    auto blocks = std::vector<::gemmi::cif::Block*>{};
    for (auto& block : document.blocks) {
        if (block.has_tag("_atom_site.id")) {
            blocks.push_back(std::addressof(block));
        }
    }
    if (every_block && blocks.size() != inputs.size()) {
        throw OutputError{OutputErrorKind::target_mismatch,
                          "Gemmi target molecule count does not match the calculation input"};
    }
    auto charged = std::vector<bool>(blocks.size(), false);
    for (std::size_t index = 0; index < inputs.size(); ++index) {
        const auto block_index = inputs[index].identity.record_index;
        if (block_index >= blocks.size() || charged[block_index]) {
            throw OutputError{OutputErrorKind::target_mismatch,
                              "Gemmi target molecules do not match the calculation input", index};
        }
        charged[block_index] = true;
        auto& block = *blocks[block_index];
        try {
            if (existing == ExistingCharges::error &&
                (block.has_mmcif_category(metadata_category) ||
                 block.has_mmcif_category(charges_category))) {
                throw OutputError{OutputErrorKind::existing_charges,
                                  "Gemmi target already contains partial charge categories"};
            }
            const auto mapping = attached_mapping(block, inputs[index]);
            append_charge_sets(block, results, index, mapping.atom_site_ids,
                               existing == ExistingCharges::append);
        } catch (...) {
            rethrow_for_molecule(index);
        }
    }
    ::gemmi::cif::write_cif_to_stream(output, document);
    if (!output) {
        throw std::runtime_error{"failed to write mmCIF output"};
    }
}

} // namespace

auto existing_charges_from_string(const std::string_view value) -> ExistingCharges {
    if (value == "error") {
        return ExistingCharges::error;
    }
    if (value == "replace") {
        return ExistingCharges::replace;
    }
    if (value == "append") {
        return ExistingCharges::append;
    }
    throw std::invalid_argument{"unknown policy for existing charges: " + std::string{value}};
}

MmcifWriter::MmcifWriter(std::ostream& output) : output_{std::addressof(output)} {}

auto MmcifWriter::write(const ChargeCalculationResult& result) const -> void {
    const auto* const single = std::addressof(result);
    write(ResultSpan{&single, 1});
}

auto MmcifWriter::write(const ResultSpan results) const -> void {
    validate_results_output(results);
    auto document = ::gemmi::cif::Document{};
    const auto inputs = results.front()->inputs();
    auto names = BlockNames{};
    document.blocks.reserve(inputs.size());
    for (std::size_t index = 0; index < inputs.size(); ++index) {
        // Names are already unique, so the document's own duplicate search is not needed.
        auto& block = document.blocks.emplace_back(names.unique(block_name(inputs[index], index)));
        try {
            write_result_block(block, results, index);
        } catch (...) {
            rethrow_for_molecule(index);
        }
    }
    ::gemmi::cif::write_cif_to_stream(*output_, document);
    if (!*output_) {
        throw std::runtime_error{"failed to write mmCIF output"};
    }
}

auto MmcifWriter::annotate(const ResultSpan results, const ExistingCharges existing) const -> void {
    // The source first: input that is not one structure file is never annotated, whatever else
    // the format could not hold of it.
    require_results(results);
    // Equal records can come from inputs that differ outside the selected sites, such as in an
    // assembly, so every record of every result must have the same source text.
    const auto first_inputs = results.front()->inputs();
    const auto& first_metadata = first_inputs.front().import_metadata;
    const auto text = first_metadata.has_value() ? first_metadata->source_text : nullptr;
    const auto same_source = [&text](const ImportedMoleculeRecord& record) {
        if (!record.import_metadata.has_value()) {
            return false;
        }
        const auto& other = record.import_metadata->source_text;
        return other == text || (other != nullptr && *other == *text);
    };
    if (text == nullptr || !std::ranges::all_of(results, [&](const auto* result) {
            return std::ranges::all_of(result->inputs(), same_source);
        })) {
        throw OutputError{OutputErrorKind::unsupported_source,
                          "mmCIF annotation requires molecules read from one PDB or mmCIF input"};
    }
    validate_results_output(results);
    if (first_metadata->format == MolecularSourceFormat::pdb) {
        // A PDB input is one molecule and has no charges; the document is its conversion.
        if (first_inputs.size() != 1) {
            throw OutputError{OutputErrorKind::unsupported_shape,
                              "mmCIF annotation of PDB input requires results for one molecule, "
                              "but the input has " +
                                  std::to_string(first_inputs.size())};
        }
        auto converted = gemmi::pdb_conversion::convert(*text, first_inputs.front());
        append_charge_sets(converted.document.blocks.front(), results, 0, converted.atom_site_ids);
        ::gemmi::cif::write_cif_to_stream(*output_, converted.document);
        if (!*output_) {
            throw std::runtime_error{"failed to write mmCIF output"};
        }
        return;
    }
    write_charged(*output_, results,
                  ::gemmi::cif::read_memory(text->data(), text->size(), "source"), existing, false);
}

// Declared in mmcif_attachment.h, outside the public namespace. The target is the caller's
// document, so a coordinate block without a calculated molecule means that it has changed.
auto gemmi::mmcif_output::write_attached(std::ostream& output, const ResultSpan results,
                                         ::gemmi::cif::Document document,
                                         const ExistingCharges existing) -> void {
    validate_results_output(results);
    write_charged(output, results, std::move(document), existing, true);
}

} // namespace chargefw::io
