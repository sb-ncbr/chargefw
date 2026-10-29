#include <chargefw/adapters/gemmi/mmcif_output.h>

#include "adapters/gemmi/mmcif_attachment.h"
#include "adapters/gemmi/mmcif_labels.h"

#include <chargefw/core/periodic_table.h>

#include <gemmi/cif.hpp>
#include <gemmi/numb.hpp>
#include <gemmi/to_cif.hpp>

#include <algorithm>
#include <cctype>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <format>
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
#include <unordered_set>
#include <utility>
#include <vector>

namespace chargefw::adapters::gemmi::mmcif_output {
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
    auto candidate = record.identity.record_id.empty() ? std::string{record.molecule.name()}
                                                       : record.identity.record_id.display_string();
    if (candidate.empty()) {
        candidate = "molecule_" + std::to_string(index + 1);
    }
    for (auto& character : candidate) {
        if (std::isalnum(static_cast<unsigned char>(character)) == 0 && character != '_' &&
            character != '-' && character != '.') {
            character = '_';
        }
    }
    return candidate;
}

[[nodiscard]] auto unique_block_name(const ::gemmi::cif::Document& document, std::string candidate)
    -> std::string {
    if (document.find_block(candidate) == nullptr) {
        return candidate;
    }
    const auto base = candidate;
    for (std::size_t suffix = 2;; ++suffix) {
        candidate = base + "_" + std::to_string(suffix);
        if (document.find_block(candidate) == nullptr) {
            return candidate;
        }
    }
}

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

[[nodiscard]] auto known_label(const std::optional<std::string>& value)
    -> std::optional<std::string_view> {
    if (!value.has_value() || value->empty() || *value == "." || *value == "?") {
        return std::nullopt;
    }
    return *value;
}

[[nodiscard]] auto label_or(const std::optional<std::string>& preferred,
                            const std::optional<std::string>& alternate,
                            const std::string_view fallback) -> std::string {
    if (const auto value = known_label(preferred); value.has_value()) {
        return std::string{*value};
    }
    if (const auto value = known_label(alternate); value.has_value()) {
        return std::string{*value};
    }
    return std::string{fallback};
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
    for (const auto& assignment : result.execution().charges->assignments()) {
        if (assignment.target.molecule_index != molecule_index) {
            continue;
        }
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

[[nodiscard]] auto attached_mapping(::gemmi::cif::Block& block,
                                    const ImportedMoleculeRecord& record) -> BlockMapping {
    if (!record.import_metadata.has_value() ||
        record.import_metadata->format != MolecularSourceFormat::mmcif) {
        throw std::invalid_argument{
            "charge attachment requires unchanged mmCIF-imported calculation inputs"};
    }
    if (record.identity.record_id.empty() ||
        record.identity.record_id.display_string() != block.name) {
        throw std::invalid_argument{
            "Gemmi target block identity does not match the calculation input"};
    }
    auto atom_sites = mmcif_labels::source_atom_sites(block);
    auto type_symbols = block.find("_atom_site.", {"type_symbol"});
    auto identity_sites =
        block.find("_atom_site.", {"Cartn_x", "Cartn_y", "Cartn_z", "?pdbx_formal_charge"});
    if (atom_sites.length() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::invalid_argument{"Gemmi target contains too many atom sites"};
    }
    if (type_symbols.length() != atom_sites.length() ||
        identity_sites.length() != atom_sites.length()) {
        throw std::invalid_argument{
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
                throw std::invalid_argument{
                    "Gemmi target site mapping does not match the calculation input"};
            }
            const auto row_index = static_cast<int>(site.position);
            auto target_id = ::gemmi::cif::as_string(atom_sites[row_index][0]);
            if (target_id != *site.id || !mapped_ids.insert(target_id).second) {
                throw std::invalid_argument{
                    "Gemmi target site mapping does not match the calculation input"};
            }
            if (!canonical_charge_site_id(target_id)) {
                throw std::invalid_argument{
                    "Gemmi target atom IDs are not representable by the charge dictionary"};
            }
            const auto& target_site = atom_sites[row_index];
            if (::gemmi::cif::as_string(type_symbols[row_index][0]) !=
                core::element_symbol(record.molecule.atom(atom_index).atomic_number())) {
                throw std::invalid_argument{
                    "Gemmi target site mapping does not match the calculation input"};
            }
            if (!site.structural_labels.has_value() ||
                mmcif_labels::decode_source_labels(target_site) != *site.structural_labels) {
                throw std::invalid_argument{
                    "Gemmi target site mapping does not match the calculation input"};
            }
            if (target_site.has(mmcif_labels::model_id_column) &&
                ::gemmi::cif::as_string(target_site[mmcif_labels::model_id_column]) != model_id) {
                throw std::invalid_argument{
                    "Gemmi target model mapping does not match the calculation input"};
            }
            const auto& position = record.molecule.conformer(conformer_index)[atom_index];
            const auto& identity_site = identity_sites[row_index];
            const auto formal_charge =
                identity_site.has(3) ? ::gemmi::cif::as_int(identity_site[3], 0) : 0;
            if (::gemmi::cif::as_number(identity_site[0]) != position.x ||
                ::gemmi::cif::as_number(identity_site[1]) != position.y ||
                ::gemmi::cif::as_number(identity_site[2]) != position.z ||
                formal_charge != record.molecule.atom(atom_index).formal_charge()) {
                throw std::invalid_argument{
                    "Gemmi target site mapping does not match the calculation input"};
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
    return reference.structural_labels.has_value() ? std::addressof(*reference.structural_labels)
                                                   : nullptr;
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
        return std::string{known_label(value).value_or(fallback)};
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

// The residue as users know it from PDB files: name, author chain, and author number.
[[nodiscard]] auto residue_description(const OutputSiteLabels& labels) -> std::string {
    return std::format("residue {} {} {}{}", labels.label_residue, labels.author_chain,
                       labels.author_sequence,
                       labels.insertion_code == "." ? "" : labels.insertion_code);
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
    std::set<std::string> atom_names;
    std::set<AtomNamePair> bonds;
};

struct ComponentBond {
    std::string component;
    std::string first_atom;
    std::string second_atom;
    core::BondOrder order;
};

struct RecordBonds {
    std::vector<ComponentBond> component_bonds;
    std::vector<const core::Bond*> connections;
};

// Splits the bonds into component bonds, written per residue name, and connections between
// residues. A component bond applies to every residue with that name, so all such residues must
// have the same bonds between the atoms they share.
[[nodiscard]] auto record_bonds(const core::Molecule& molecule,
                                const std::vector<OutputSiteLabels>& labels,
                                const std::vector<std::size_t>& residue_of) -> RecordBonds {
    auto result = RecordBonds{};
    auto residues = std::vector<ResidueBonds>(std::ranges::max(residue_of) + 1);
    for (std::size_t atom_index = 0; atom_index < molecule.atom_count(); ++atom_index) {
        const auto& site = labels[atom_index];
        auto& residue = residues[residue_of[atom_index]];
        if (residue.atom_names.empty()) {
            residue.component = site.label_residue;
            residue.description = residue_description(site);
        }
        if (!residue.atom_names.insert(site.label_atom).second) {
            throw std::invalid_argument{std::format(
                "mmCIF bond output requires unique atom names, but {} has atom '{}' more than once",
                residue.description, site.label_atom)};
        }
    }

    auto definitions = std::map<std::string, std::map<AtomNamePair, core::BondOrder>>{};
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
        residue.bonds.insert(key);
        const auto [definition, inserted] =
            definitions[residue.component].try_emplace(key, bond.order());
        if (inserted) {
            result.component_bonds.push_back({.component = residue.component,
                                              .first_atom = first,
                                              .second_atom = second,
                                              .order = bond.order()});
        } else if (definition->second != bond.order()) {
            throw std::invalid_argument{std::format(
                "mmCIF bond output cannot write bond {}-{} of {}: another {} residue has a "
                "different bond order",
                first, second, residue.description, residue.component)};
        }
    }

    // A residue without a bond that its component defines would gain that bond in readers.
    for (const auto& residue : residues) {
        for (const auto& [key, order] : definitions[residue.component]) {
            if (residue.atom_names.contains(key.first) && residue.atom_names.contains(key.second) &&
                !residue.bonds.contains(key)) {
                throw std::invalid_argument{std::format(
                    "mmCIF bond output cannot write bond {}-{} of {}: {} has both atoms but no "
                    "bond between them",
                    key.first, key.second, residue.component, residue.description)};
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
                         std::string{mmcif_bond_order(bond.order)}, std::to_string(ordinal++)});
    }
}

// Appends the _struct_conn partner columns of one atom, in the order of connection_columns().
auto append_connection_partner(std::vector<std::string>& row, const OutputSiteLabels& site)
    -> void {
    // Readers find a partner by author chain and number, or by label chain and number.
    if (site.author_sequence == "." && site.label_sequence == ".") {
        throw std::invalid_argument{"mmCIF bond output cannot write a bond to " +
                                    residue_description(site) + ": it has no residue number"};
    }
    // Readers such as Gemmi name atoms by their author names but look up partners by the label
    // names, so a partner whose names differ would not be found.
    if (site.label_atom != site.author_atom || site.label_residue != site.author_residue) {
        throw std::invalid_argument{std::format(
            "mmCIF bond output cannot write a bond to atom {} of {}: its author and label atom or "
            "residue names differ",
            site.label_atom, residue_description(site))};
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
        append_connection_partner(row, labels[bond->first_atom_index()]);
        append_connection_partner(row, labels[bond->second_atom_index()]);
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

auto validate_result_output(const ChargeCalculationResult& result) -> void {
    if (!result.execution().calculated() || !result.execution().charges.has_value()) {
        throw std::invalid_argument{"mmCIF output requires a successful calculation"};
    }
    if (result.inputs().empty()) {
        throw std::invalid_argument{"mmCIF output requires at least one molecule record"};
    }
    for (const auto& record : result.inputs()) {
        const auto& molecule = record.molecule;
        if (molecule.atom_count() == 0) {
            throw std::invalid_argument{"mmCIF output requires at least one atom per record"};
        }
        if (molecule.conformer_count() == 0) {
            throw std::invalid_argument{"mmCIF output requires coordinates"};
        }
        if (molecule.atom_count() != 0 &&
            molecule.conformer_count() >
                static_cast<std::size_t>(std::numeric_limits<int>::max()) / molecule.atom_count()) {
            throw std::invalid_argument{"mmCIF output contains too many atom sites"};
        }
        for (const auto& conformer : molecule.conformers()) {
            for (const auto& position : conformer.positions()) {
                if (!std::isfinite(position.x) || !std::isfinite(position.y) ||
                    !std::isfinite(position.z)) {
                    throw std::invalid_argument{"mmCIF output coordinates must be finite"};
                }
            }
        }
    }
    for (const auto& assignment : result.execution().charges->assignments()) {
        for (const auto charge : assignment.charges.values()) {
            if (charge < -5.0 || charge > 5.0) {
                throw std::invalid_argument{"mmCIF charge is outside the dictionary range [-5, 5]"};
            }
        }
    }
}

using ResultSpan = std::span<const ChargeCalculationResult* const>;

// Several charge sets share one file only when they were calculated for the same molecules and
// can be told apart; the file records method and parameter set but not options or execution.
auto validate_results_output(const ResultSpan results) -> void {
    if (results.empty()) {
        throw std::invalid_argument{"mmCIF output requires at least one calculation result"};
    }
    auto charge_sets = std::set<std::pair<std::string, std::optional<std::string>>>{};
    for (const auto* result : results) {
        if (result == nullptr) {
            throw std::invalid_argument{"mmCIF output requires calculation results, not null"};
        }
        validate_result_output(*result);
        if (result != results.front() &&
            !std::ranges::equal(result->inputs(), results.front()->inputs())) {
            throw std::invalid_argument{
                "mmCIF output requires charge sets calculated for the same molecules"};
        }
        const auto& charge_set = *result->execution().charges;
        const auto parameter_set = charge_set.parameter_set_id().transform(
            [](const std::string_view value) { return std::string{value}; });
        if (!charge_sets.emplace(std::string{charge_set.method_id()}, parameter_set).second) {
            throw std::invalid_argument{std::format(
                "mmCIF output cannot hold two charge sets of method {} with parameter set {}: "
                "the file cannot distinguish their options or execution",
                charge_set.method_id(), parameter_set.value_or("none"))};
        }
    }
}

auto append_charge_metadata(::gemmi::cif::Table& metadata, const std::string_view id,
                            const charges::ChargeSet& charge_set,
                            const std::string_view generator_version) -> void {
    const auto parameter_set = charge_set.parameter_set_id();
    metadata.append_row({std::string{id}, "empirical", quote(charge_set.method_id()),
                         parameter_set.has_value() ? quote(*parameter_set) : ".", quote("ChargeFW"),
                         quote(generator_version.empty() ? "unknown" : generator_version)});
}

// Appends one metadata row per charge set and its charges for one record. Charge set k has id
// k + 1 and covers every conformer; site_ids[conformer][atom] are the record's _atom_site.id
// values.
auto append_charge_sets(::gemmi::cif::Block& block, const ResultSpan results,
                        const std::size_t record_index,
                        const std::vector<std::vector<std::string>>& site_ids,
                        const std::string_view generator_version) -> void {
    ensure_dictionary(block);
    block.init_mmcif_loop(metadata_category, metadata_columns());
    block.init_mmcif_loop(charges_category, charge_columns());

    // Category initialization can invalidate tables, so acquire them only after both mutations.
    auto metadata = block.find(metadata_category, metadata_columns());
    auto charge_rows = block.find(charges_category, charge_columns());
    for (std::size_t result_index = 0; result_index < results.size(); ++result_index) {
        const auto& result = *results[result_index];
        const auto id = std::to_string(result_index + 1);
        append_charge_metadata(metadata, id, *result.execution().charges, generator_version);
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
                        const std::size_t record_index, const std::string_view generator_version)
    -> void {
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
    append_charge_sets(block, results, record_index, site_ids, generator_version);
    append_bonds(block, record, generated_atom_ids);
}

// Declares the charges dictionary once. Declarations of other charges dictionary versions, such as
// mmcif_charges_v10.dic in older files, no longer describe the charges and are removed; other
// dictionaries stay.
auto ensure_dictionary(::gemmi::cif::Block& block) -> void {
    auto table =
        block.find_or_add("_audit_conform.", {"dict_name", "dict_version", "dict_location"});
    table.ensure_loop();
    for (auto index = static_cast<int>(table.length()); index-- > 0;) {
        if (::gemmi::cif::as_string(table[index][0]).starts_with("mmcif_charges_")) {
            table.remove_row(index);
        }
    }
    table.append_row({dictionary_name, dictionary_version, dictionary_location});
}

} // namespace

MmcifWriter::MmcifWriter(std::ostream& output) : output_{std::addressof(output)} {}

auto MmcifWriter::write(const ChargeCalculationResult& result,
                        const std::string_view generator_version) const -> void {
    const auto* const single = std::addressof(result);
    write(ResultSpan{&single, 1}, generator_version);
}

auto MmcifWriter::write(const ResultSpan results, const std::string_view generator_version) const
    -> void {
    validate_results_output(results);
    auto document = ::gemmi::cif::Document{};
    const auto inputs = results.front()->inputs();
    for (std::size_t index = 0; index < inputs.size(); ++index) {
        auto& block =
            document.add_new_block(unique_block_name(document, block_name(inputs[index], index)));
        write_result_block(block, results, index, generator_version);
    }
    ::gemmi::cif::write_cif_to_stream(*output_, document);
    if (!*output_) {
        throw std::runtime_error{"failed to write mmCIF output"};
    }
}

auto write_attached(std::ostream& output, const ResultSpan results,
                    const ::gemmi::cif::Document& source, const bool overwrite,
                    const std::string_view generator_version) -> void {
    validate_results_output(results);
    const auto inputs = results.front()->inputs();
    auto document = source;
    auto blocks = std::vector<::gemmi::cif::Block*>{};
    for (auto& block : document.blocks) {
        if (block.has_tag("_atom_site.id")) {
            blocks.push_back(std::addressof(block));
        }
    }
    if (blocks.size() != inputs.size()) {
        throw std::invalid_argument{
            "Gemmi target molecule count does not match the calculation input"};
    }
    for (std::size_t index = 0; index < blocks.size(); ++index) {
        auto& block = *blocks[index];
        if (!overwrite && (block.has_mmcif_category(metadata_category) ||
                           block.has_mmcif_category(charges_category))) {
            throw std::invalid_argument{"Gemmi target already contains partial charge categories"};
        }
        const auto mapping = attached_mapping(block, inputs[index]);
        append_charge_sets(block, results, index, mapping.atom_site_ids, generator_version);
    }
    ::gemmi::cif::write_cif_to_stream(output, document);
    if (!output) {
        throw std::runtime_error{"failed to write mmCIF output"};
    }
}

} // namespace chargefw::adapters::gemmi::mmcif_output
