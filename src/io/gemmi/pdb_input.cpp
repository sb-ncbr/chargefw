#include <chargefw/io/pdb_input.h>

#include "bonds.h"
#include "io/text_encoding.h"
#include "pdb_conversion.h"
#include "selection.h"
#include "structure_import.h"

#include <chargefw/core/bond.h>

#include <gemmi/align.hpp>
#include <gemmi/enumstr.hpp>
#include <gemmi/pdb.hpp>
#include <gemmi/polyheur.hpp>
#include <gemmi/to_mmcif.hpp>

#include <algorithm>
#include <array>
#include <cctype>
#include <charconv>
#include <cstddef>
#include <cstdlib>
#include <format>
#include <functional>
#include <istream>
#include <iterator>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace chargefw::io {
namespace {

[[nodiscard]] auto optional_text(const std::string_view value) -> std::optional<std::string> {
    return value.empty() ? std::nullopt : std::optional<std::string>{value};
}

[[nodiscard]] auto optional_character(const char value) -> std::optional<std::string> {
    return value == '\0' || value == ' ' ? std::nullopt
                                         : std::optional<std::string>{std::string{value}};
}

[[nodiscard]] auto field(const std::string_view line, const std::size_t offset,
                         const std::size_t width) -> std::optional<std::string> {
    if (offset >= line.size()) {
        return std::nullopt;
    }
    const auto value = line.substr(offset, std::min(width, line.size() - offset));
    const auto first = value.find_first_not_of(' ');
    if (first == std::string_view::npos) {
        return std::nullopt;
    }
    const auto last = value.find_last_not_of(' ');
    return std::string{value.substr(first, last - first + 1)};
}

[[nodiscard]] auto normalized_sequence(const std::optional<std::string>& sequence)
    -> std::optional<std::string> {
    if (!sequence.has_value()) {
        return std::nullopt;
    }
    if (sequence->front() < 'A') {
        auto value = int{};
        const auto [end, error] =
            std::from_chars(sequence->data(), sequence->data() + sequence->size(), value);
        if (error != std::errc{} || end != sequence->data() + sequence->size()) {
            throw std::runtime_error{"PDB source contains an invalid residue sequence ID"};
        }
        return std::to_string(value);
    }

    const auto value = std::strtol(sequence->c_str(), nullptr, 36) - 466560L + 10000L;
    return std::to_string(value);
}

struct PdbSourceSite {
    SourceAtomReference reference;
    SourceStructuralLabels identity;
};

struct PdbSourceModel {
    std::optional<std::string> id;
    std::vector<PdbSourceSite> sites;
};

struct PdbSource {
    std::vector<PdbSourceModel> models;
    SourceConnectivity connectivity = SourceConnectivity::absent;
};

struct StructuralLabelsLess {
    [[nodiscard]] auto operator()(const SourceStructuralLabels& lhs,
                                  const SourceStructuralLabels& rhs) const -> bool {
        return std::tie(lhs.author.atom, lhs.author.residue, lhs.author.chain, lhs.author.sequence,
                        lhs.label.atom, lhs.label.residue, lhs.label.chain, lhs.label.sequence,
                        lhs.entity, lhs.insertion_code, lhs.alternate_location, lhs.segment) <
               std::tie(rhs.author.atom, rhs.author.residue, rhs.author.chain, rhs.author.sequence,
                        rhs.label.atom, rhs.label.residue, rhs.label.chain, rhs.label.sequence,
                        rhs.entity, rhs.insertion_code, rhs.alternate_location, rhs.segment);
    }
};

// Reads model IDs and connectivity from the whole file, but atom sites only for the first
// retained_models models; the others are not mapped.
[[nodiscard]] auto parse_source(const std::string& contents, const std::size_t retained_models)
    -> PdbSource {
    auto result = PdbSource{};
    auto* current = static_cast<PdbSourceModel*>(nullptr);
    auto source_position = std::size_t{0};
    auto input = std::istringstream{contents};
    auto line = std::string{};
    while (std::getline(input, line)) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        const auto view = std::string_view{line};
        auto record = std::string{view.substr(0, std::min(view.size(), std::size_t{6}))};
        std::ranges::transform(record, record.begin(), [](const unsigned char character) {
            return static_cast<char>(std::toupper(character));
        });
        if (record.starts_with("MODEL ")) {
            result.models.push_back(PdbSourceModel{.id = field(view, 10, 4), .sites = {}});
            current = std::addressof(result.models.back());
            continue;
        }
        if (record.starts_with("ENDMDL")) {
            current = nullptr;
            continue;
        }
        if (record.starts_with("END")) {
            break;
        }
        if (record.starts_with("CONECT") || record.starts_with("LINK  ") ||
            record.starts_with("SSBOND")) {
            result.connectivity = SourceConnectivity::present;
        }
        if (!record.starts_with("ATOM  ") && !record.starts_with("HETATM")) {
            continue;
        }
        if (current == nullptr) {
            result.models.push_back(PdbSourceModel{.id = std::nullopt, .sites = {}});
            current = std::addressof(result.models.back());
        }
        if (result.models.size() > retained_models) {
            continue;
        }
        auto labels =
            SourceStructuralLabels{.author = SourceHierarchyLabels{.atom = field(view, 12, 4),
                                                                   .residue = field(view, 17, 3),
                                                                   .chain = field(view, 20, 2),
                                                                   .sequence = field(view, 22, 4)},
                                   .label = {},
                                   .entity = std::nullopt,
                                   .insertion_code = field(view, 26, 1),
                                   .alternate_location = field(view, 16, 1),
                                   .segment = field(view, 72, 4)};
        auto identity = labels;
        identity.author.sequence = normalized_sequence(labels.author.sequence);
        current->sites.push_back(PdbSourceSite{
            .reference = SourceAtomReference{.position = source_position++,
                                             .id = field(view, 6, 5),
                                             .structural_labels =
                                                 std::make_shared<const SourceStructuralLabels>(
                                                     std::move(labels))},
            .identity = std::move(identity)});
    }
    return result;
}

[[nodiscard]] auto selected_labels(const ::gemmi::const_CRA& site) -> SourceStructuralLabels {
    return SourceStructuralLabels{
        .author = SourceHierarchyLabels{.atom = optional_text(site.atom->name),
                                        .residue = optional_text(site.residue->name),
                                        .chain = optional_text(site.chain->name),
                                        .sequence = optional_text(site.residue->seqid.num.str())},
        .label = {},
        .entity = std::nullopt,
        .insertion_code = optional_character(site.residue->seqid.icode),
        .alternate_location = optional_character(site.atom->altloc),
        .segment = optional_text(site.residue->segment)};
}

[[nodiscard]] auto
source_models(const ::gemmi::Structure& structure, const StructureOptions options,
              const std::span<const gemmi::selection::SelectedModel> selected_models,
              const std::vector<PdbSourceModel>& source)
    -> std::vector<gemmi::structure_import::SourceModelMapping> {
    const auto retained_count =
        options.conformers == ConformerSelection::all ? structure.models.size() : std::size_t{1};
    if (source.size() != structure.models.size()) {
        throw std::runtime_error{"PDB source mapping does not match parsed models"};
    }
    auto result = std::vector<gemmi::structure_import::SourceModelMapping>{};
    result.reserve(retained_count);
    for (std::size_t model_index = 0; model_index < retained_count; ++model_index) {
        const auto& selected = selected_models[model_index];
        auto sites = std::vector<SourceAtomReference>{};
        sites.reserve(selected.sites().size());
        auto source_indices =
            std::map<SourceStructuralLabels, std::vector<std::size_t>, StructuralLabelsLess>{};
        for (std::size_t index = 0; index < source[model_index].sites.size(); ++index) {
            source_indices[source[model_index].sites[index].identity].push_back(index);
        }
        for (const auto& site : selected.sites()) {
            const auto labels = selected_labels(site);
            const auto found = source_indices.find(labels);
            if (found == source_indices.end()) {
                throw std::runtime_error{"PDB source mapping does not match selected atoms"};
            }
            if (found->second.size() != 1) {
                throw std::runtime_error{"PDB source contains ambiguous selected atom identity"};
            }
            sites.push_back(source[model_index].sites[found->second.front()].reference);
            source_indices.erase(found);
        }
        result.push_back(gemmi::structure_import::SourceModelMapping{
            .conformer = SourceConformerReference{
                .position = model_index, .id = source[model_index].id, .sites = std::move(sites)}});
    }
    return result;
}

} // namespace

PdbReader::PdbReader(std::istream& input, std::string source_name,
                     const ::chargefw::io::StructureOptions options) {
    std::string contents{std::istreambuf_iterator<char>{input}, std::istreambuf_iterator<char>{}};
    if (input.bad()) {
        throw std::runtime_error{"failed to read PDB input"};
    }
    if (const auto line = io::detail::first_non_ascii_line(contents)) {
        throw std::runtime_error{"line " + std::to_string(*line) +
                                 " contains a non-ASCII character; PDB files must be ASCII text"};
    }

    auto structure = ::gemmi::read_pdb_string(contents, source_name);
    if (structure.models.empty()) {
        throw std::runtime_error{"structural input contains no models"};
    }
    ::gemmi::add_entity_types(structure, false);
    const auto name = structure.name.empty() ? source_name : structure.name;
    auto selected_models =
        gemmi::selection::select_models(structure, options.selection, options.conformers);
    const auto retained_models =
        options.conformers == ConformerSelection::all ? structure.models.size() : std::size_t{1};
    const auto parsed_source = parse_source(contents, retained_models);
    auto mappings = source_models(structure, options, selected_models, parsed_source.models);
    for (std::size_t index = 0; index < selected_models.size(); ++index) {
        selected_models[index].normalize_component_names(mappings[index].conformer.sites);
    }
    auto explicit_bonds = gemmi::bonds::ExplicitBonds{};
    if (options.bonds == BondStrategy::explicit_bonds || options.bonds == BondStrategy::hybrid) {
        explicit_bonds = gemmi::bonds::explicit_pdb(structure, selected_models.front());
    }
    record_ = gemmi::structure_import::make_record(
        structure, selected_models,
        MoleculeRecordIdentity{
            .source_name = std::move(source_name), .record_index = 0, .record_id = name},
        options.selection, options.bonds, options.conformers, std::move(explicit_bonds.bonds), name,
        MolecularSourceFormat::pdb, std::move(mappings), parsed_source.connectivity);
    record_->import_metadata->source_text =
        std::make_shared<const std::string>(std::move(contents));
    if (explicit_bonds.skipped_symmetry_connections > 0) {
        record_->import_warnings.push_back(
            gemmi::bonds::symmetry_connections_issue(explicit_bonds.skipped_symmetry_connections));
    }
}

auto PdbReader::next() -> std::optional<ImportedMoleculeRecord> {
    return std::exchange(record_, std::nullopt);
}

namespace {

// A site of the structure read again, in the order in which Gemmi writes atom sites.
struct ConvertedSite {
    ::gemmi::const_CRA site;
    std::size_t number = 0; // Counted through all models; the _atom_site.id is number + 1.
};

[[nodiscard]] auto model_sites(const ::gemmi::Model& model, std::size_t number)
    -> std::vector<ConvertedSite> {
    auto result = std::vector<ConvertedSite>{};
    for (const auto& chain : model.chains) {
        for (const auto& residue : chain.residues) {
            for (const auto& atom : residue.atoms) {
                result.push_back({.site = {.chain = std::addressof(chain),
                                           .residue = std::addressof(residue),
                                           .atom = std::addressof(atom)},
                                  .number = number++});
            }
        }
    }
    return result;
}

// Finds each calculated atom of one model by the labels its reader matched it by, and checks its
// element and coordinates.
[[nodiscard]] auto calculated_sites(const std::span<const ConvertedSite> sites,
                                    const SourceConformerReference& conformer,
                                    const core::Molecule& molecule,
                                    const std::size_t conformer_index)
    -> std::vector<const ConvertedSite*> {
    auto by_labels =
        std::map<SourceStructuralLabels, std::vector<const ConvertedSite*>, StructuralLabelsLess>{};
    for (const auto& site : sites) {
        by_labels[selected_labels(site.site)].push_back(std::addressof(site));
    }
    auto result = std::vector<const ConvertedSite*>{};
    result.reserve(conformer.sites.size());
    for (std::size_t index = 0; index < conformer.sites.size(); ++index) {
        const auto& labels = conformer.sites[index].structural_labels;
        if (labels == nullptr) {
            throw std::invalid_argument{"PDB annotation requires the labels of the PDB input"};
        }
        auto identity = *labels;
        identity.author.sequence = normalized_sequence(labels->author.sequence);
        const auto found = by_labels.find(identity);
        if (found == by_labels.end() || found->second.size() != 1) {
            throw std::invalid_argument{"PDB source does not match the calculation input"};
        }
        const auto& atom = *found->second.front()->site.atom;
        const auto& position = molecule.conformer(conformer_index)[index];
        if (atom.element.atomic_number() != molecule.atom(index).atomic_number() ||
            atom.pos.x != position.x || atom.pos.y != position.y || atom.pos.z != position.z) {
            throw std::invalid_argument{"PDB source does not match the calculation input"};
        }
        result.push_back(found->second.front());
    }
    return result;
}

// A connection of the first model between two of its sites, as it is written to _struct_conn.
struct OutputConnection {
    std::string id;
    ::gemmi::Connection::Type type = ::gemmi::Connection::Covale;
    const ConvertedSite* first = nullptr;
    const ConvertedSite* second = nullptr;
    // Empty when the PDB file does not say; the readers then take the partners as local.
    std::string first_symmetry = {};
    std::string second_symmetry = {};
    bool to_symmetry_copy = false;
    // The distance a LINK or SSBOND record reports; zero when it reports none.
    double reported_distance = 0.0;

    // A bond between atoms of the file, which the readers add to the molecule.
    [[nodiscard]] auto local_bond() const -> bool {
        return !to_symmetry_copy &&
               (type == ::gemmi::Connection::Covale || type == ::gemmi::Connection::Disulf);
    }
};

using SitePair = std::pair<const ConvertedSite*, const ConvertedSite*>;

[[nodiscard]] auto site_pair(const ConvertedSite* first, const ConvertedSite* second) -> SitePair {
    return std::less<>{}(second, first) ? SitePair{second, first} : SitePair{first, second};
}

// The LINK and SSBOND connections whose atoms the PDB reader finds, with the symmetry operators
// of the file relative to the first partner: 1655 and 1755 become 1_555 and 1_655, and equal
// operators are a local connection. Different operators whose first is not a translation, such
// as 2555 and 3555, cannot be stated that way and are refused.
[[nodiscard]] auto
link_connections(const ::gemmi::Structure& structure,
                 const std::unordered_map<const ::gemmi::Atom*, const ConvertedSite*>& sites)
    -> std::vector<OutputConnection> {
    auto result = std::vector<OutputConnection>{};
    const auto partners = gemmi::bonds::ResidueIndex{structure.models.front()};
    for (const auto& connection : structure.connections) {
        const auto first = sites.find(gemmi::bonds::partner_atom(partners, connection.partner1));
        const auto second = sites.find(gemmi::bonds::partner_atom(partners, connection.partner2));
        if (first == sites.end() || second == sites.end()) {
            continue;
        }
        auto output = OutputConnection{.id = connection.name,
                                       .type = connection.type,
                                       .first = first->second,
                                       .second = second->second,
                                       .reported_distance = connection.reported_distance};
        if (connection.asu == ::gemmi::Asu::Same) {
            output.first_symmetry = output.second_symmetry = "1_555";
        } else if (connection.asu == ::gemmi::Asu::Different) {
            const auto& symmetry = connection.reported_sym;
            const auto shift = [&symmetry](const std::size_t axis) { return 5 + symmetry[axis]; };
            if (symmetry[0] <= 0 || symmetry[0] == 99 ||
                std::ranges::any_of(std::array<std::size_t, 3>{1, 2, 3}, [&](const auto axis) {
                    return shift(axis) < 0 || shift(axis) > 9;
                })) {
                throw std::invalid_argument{std::format(
                    "PDB connection {} has symmetry operators that mmCIF output cannot state",
                    connection.name)};
            }
            output.to_symmetry_copy = true;
            output.first_symmetry = "1_555";
            output.second_symmetry =
                std::format("{}_{}{}{}", symmetry[0], shift(1), shift(2), shift(3));
        }
        result.push_back(std::move(output));
    }
    return result;
}

// Adds the CONECT bonds of the first model as covalent connections, each pair once, unless a LINK
// or SSBOND already states a bond between the two atoms of the file.
void add_conect_connections(std::vector<OutputConnection>& connections,
                            const ::gemmi::Structure& structure,
                            const std::span<const ConvertedSite> sites) {
    if (structure.conect_map.empty()) {
        return;
    }
    auto by_serial = std::unordered_map<int, const ConvertedSite*>{};
    auto repeated = std::set<int>{};
    for (const auto& site : sites) {
        if (!by_serial.emplace(site.site.atom->serial, std::addressof(site)).second) {
            repeated.insert(site.site.atom->serial);
        }
    }
    auto joined = std::set<SitePair>{};
    for (const auto& connection : connections) {
        if (connection.local_bond()) {
            joined.insert(site_pair(connection.first, connection.second));
        }
    }
    auto count = std::size_t{0};
    for (const auto& [serial, partners] : structure.conect_map) {
        for (const auto partner : partners) {
            const auto first = by_serial.find(serial);
            const auto second = by_serial.find(partner);
            if (serial == partner || first == by_serial.end() || second == by_serial.end()) {
                continue;
            }
            if (repeated.contains(serial) || repeated.contains(partner)) {
                throw std::invalid_argument{std::format(
                    "PDB CONECT records refer to atom serial number {}, which several atoms have",
                    repeated.contains(serial) ? serial : partner)};
            }
            if (!joined.insert(site_pair(first->second, second->second)).second) {
                continue;
            }
            connections.push_back({.id = "conect" + std::to_string(++count),
                                   .first = first->second,
                                   .second = second->second,
                                   .first_symmetry = "1_555",
                                   .second_symmetry = "1_555"});
        }
    }
}

// mmCIF names an atom by chain, residue, atom name, and alternate location. Atoms that differ only
// in the PDB segment would be confused there, as atom sites and as connection partners.
class DistinctSites {
  public:
    explicit DistinctSites(const std::span<const ConvertedSite> sites) {
        for (const auto& site : sites) {
            ++counts_[address(site)];
        }
    }

    void require(const ConvertedSite& site) const {
        if (counts_.at(address(site)) != 1) {
            const auto& [chain, residue, atom] = site.site;
            throw std::invalid_argument{std::format(
                "mmCIF output cannot tell atom {} of residue {} {} {} from another atom with the "
                "same names in a different PDB segment",
                atom->name, residue->name, chain->name, residue->seqid.str())};
        }
    }

  private:
    using Address =
        std::tuple<std::string_view, int, char, std::string_view, std::string_view, char>;

    [[nodiscard]] static auto address(const ConvertedSite& site) -> Address {
        const auto& [chain, residue, atom] = site.site;
        return {chain->name,          residue->seqid.num.value,
                residue->seqid.icode, residue->name,
                atom->name,           atom->altloc};
    }

    std::map<Address, std::size_t> counts_;
};

// Writes _struct_conn from the sites themselves, so that each partner is the atom the PDB reader
// used. A connection between calculated atoms states the bond order of the record: a connection
// without one is a single bond in mmCIF, where the record may have the order of a residue
// template.
void write_connections(::gemmi::cif::Block& block,
                       const std::span<const OutputConnection> connections,
                       const DistinctSites& distinct,
                       const std::span<const ConvertedSite* const> calculated,
                       const core::Molecule& molecule) {
    if (connections.empty()) {
        return;
    }
    auto atom_indices = std::unordered_map<const ConvertedSite*, std::size_t>{};
    for (std::size_t index = 0; index < calculated.size(); ++index) {
        atom_indices.emplace(calculated[index], index);
    }
    auto bond_orders = std::map<std::pair<std::size_t, std::size_t>, core::BondOrder>{};
    for (const auto& bond : molecule.bonds()) {
        bond_orders.emplace(std::minmax(bond.first_atom_index(), bond.second_atom_index()),
                            bond.order());
    }
    const auto order = [&](const OutputConnection& connection) -> std::string {
        const auto first = atom_indices.find(connection.first);
        const auto second = atom_indices.find(connection.second);
        if (!connection.local_bond() || first == atom_indices.end() ||
            second == atom_indices.end()) {
            return "?";
        }
        const auto bond = bond_orders.find(std::minmax(first->second, second->second));
        if (bond == bond_orders.end()) {
            return "?";
        }
        return bond->second == core::BondOrder::SINGLE   ? "sing"
               : bond->second == core::BondOrder::DOUBLE ? "doub"
                                                         : "trip";
    };

    auto columns = std::vector<std::string>{"id", "conn_type_id"};
    for (const auto* const partner : {"ptnr1", "ptnr2"}) {
        for (const auto* const item :
             {"{}_label_asym_id", "{}_label_comp_id", "{}_label_seq_id", "{}_label_atom_id",
              "pdbx_{}_label_alt_id", "{}_auth_asym_id", "{}_auth_seq_id", "pdbx_{}_PDB_ins_code",
              "{}_symmetry"}) {
            columns.push_back(std::vformat(item, std::make_format_args(partner)));
        }
    }
    columns.emplace_back("pdbx_dist_value");
    columns.emplace_back("pdbx_value_order");
    auto& loop = block.init_mmcif_loop("_struct_conn.", columns);
    auto types = std::set<std::string>{};
    const auto text = [](const std::string& value, const char* const missing) {
        return value.empty() ? std::string{missing} : ::gemmi::cif::quote(value);
    };
    for (const auto& connection : connections) {
        auto row = std::vector<std::string>{text(connection.id, "?"),
                                            ::gemmi::connection_type_to_string(connection.type)};
        types.insert(row.back());
        for (const auto& [site, symmetry] :
             {std::pair{connection.first, std::cref(connection.first_symmetry)},
              std::pair{connection.second, std::cref(connection.second_symmetry)}}) {
            distinct.require(*site);
            const auto& [chain, residue, atom] = site->site;
            row.push_back(text(residue->subchain, "."));
            row.push_back(text(residue->name, "?"));
            row.push_back(residue->label_seq.str('.'));
            row.push_back(text(atom->name, "?"));
            row.emplace_back(1, atom->altloc_or('?'));
            row.push_back(text(chain->name, "?"));
            row.push_back(residue->seqid.num.str());
            row.emplace_back(1, residue->seqid.has_icode() ? residue->seqid.icode : '?');
            row.push_back(text(symmetry.get(), "?"));
        }
        row.push_back(connection.reported_distance > 0.0
                          ? std::format("{}", connection.reported_distance)
                          : "?");
        row.push_back(order(connection));
        loop.add_row(row);
    }
    auto& type_loop = block.init_mmcif_loop("_struct_conn_type.", {"id"});
    for (const auto& type : types) {
        type_loop.add_row({type});
    }
}

} // namespace

auto gemmi::pdb_conversion::convert(const std::string& text, const ImportedMoleculeRecord& record)
    -> ConvertedPdb {
    if (!record.import_metadata.has_value() ||
        record.import_metadata->format != MolecularSourceFormat::pdb) {
        throw std::invalid_argument{"PDB annotation requires a molecule read from PDB"};
    }
    const auto& metadata = *record.import_metadata;
    auto structure = ::gemmi::read_pdb_string(text, record.identity.source_name);
    if (structure.models.size() < metadata.conformers.size()) {
        throw std::invalid_argument{"PDB source does not match the calculation input"};
    }
    ::gemmi::setup_entities(structure);
    ::gemmi::assign_label_seq_id(structure, false);
    if (structure.name.empty()) {
        structure.name = "structure";
    }

    auto result = ConvertedPdb{};
    auto first_model_sites = std::vector<ConvertedSite>{};
    auto first_calculated = std::vector<const ConvertedSite*>{};
    auto number = std::size_t{0};
    for (std::size_t index = 0; index < structure.models.size(); ++index) {
        auto sites = model_sites(structure.models[index], number);
        number += sites.size();
        if (index >= metadata.conformers.size()) {
            continue;
        }
        const auto calculated =
            calculated_sites(sites, metadata.conformers[index], record.molecule, index);
        auto& ids = result.atom_site_ids.emplace_back();
        ids.reserve(calculated.size());
        for (const auto* site : calculated) {
            ids.push_back(std::to_string(site->number + 1));
        }
        if (index == 0) {
            // The sites keep their addresses when the vector is moved.
            first_model_sites = std::move(sites);
            first_calculated = calculated;
        }
    }
    auto site_of_atom = std::unordered_map<const ::gemmi::Atom*, const ConvertedSite*>{};
    for (const auto& site : first_model_sites) {
        site_of_atom.emplace(site.site.atom, std::addressof(site));
    }
    auto connections = link_connections(structure, site_of_atom);
    add_conect_connections(connections, structure, first_model_sites);

    // Gemmi's _struct_conn finds partners and symmetry images anew, which can differ from what was
    // calculated, so the connections are written here. Its _chem_comp has only the residue names
    // and, for names longer than three characters, their short forms; without those it says
    // nothing the atom sites do not.
    auto groups = ::gemmi::MmcifOutputGroups{true};
    groups.conn = false;
    groups.chem_comp = !structure.shortened_ccd_codes.empty();
    result.document = ::gemmi::make_mmcif_document(structure, groups);
    auto& block = result.document.blocks.front();
    for (auto& character : block.name) {
        if (std::isgraph(static_cast<unsigned char>(character)) == 0) {
            character = '_';
        }
    }
    const auto distinct = DistinctSites{first_model_sites};
    for (const auto* site : first_calculated) {
        distinct.require(*site);
    }
    write_connections(block, connections, distinct, first_calculated, record.molecule);
    return result;
}

} // namespace chargefw::io
