#include <chargefw/io/parse_error.h>
#include <chargefw/io/pdb_input.h>

#include "bonds.h"
#include "io/text_encoding.h"
#include "pdb_conversion.h"
#include "selection.h"
#include "structure_import.h"

#include <chargefw/core/bond.h>

#include <gemmi/align.hpp>
#include <gemmi/pdb.hpp>
#include <gemmi/polyheur.hpp>
#include <gemmi/to_mmcif.hpp>

#include <algorithm>
#include <cctype>
#include <charconv>
#include <cstddef>
#include <cstdlib>
#include <format>
#include <istream>
#include <iterator>
#include <map>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
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
        throw ParseError{"non-ASCII character; PDB files must be ASCII text",
                         {.source_name = source_name, .line = line}};
    }

    auto structure = ::gemmi::read_pdb_string(contents, source_name);
    if (structure.models.empty()) {
        throw std::runtime_error{"structural input contains no models"};
    }
    ::gemmi::add_entity_types(structure, false);
    // The ID code of the HEADER record, such as 1CBN; otherwise the file name without .pdb.
    const auto entry_id = structure.info.find("_entry.id");
    const auto name = entry_id != structure.info.end() ? entry_id->second
                      : structure.name.empty()         ? source_name
                                                       : structure.name;
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
        MolecularSourceFormat::pdb, std::move(mappings), parsed_source.connectivity,
        std::move(explicit_bonds.unknown_order_indices));
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

// The _atom_site.id of each calculated atom of one model of the structure read again. Gemmi
// writes the sites in the order of its chains, residues, and atoms, which is not always the order
// of the lines, and numbers them from 1 through all models: next_id is that of the model's first
// site and is advanced past its last. Each atom is found by the labels the reader matched it by
// and checked by its element and coordinates.
[[nodiscard]] auto calculated_site_ids(const ::gemmi::Model& model, std::size_t& next_id,
                                       const SourceConformerReference& conformer,
                                       const core::Molecule& molecule,
                                       const std::size_t conformer_index)
    -> std::vector<std::string> {
    struct Site {
        const ::gemmi::Atom* atom = nullptr;
        std::size_t id = 0;
        bool repeated = false;
    };
    auto by_labels = std::map<SourceStructuralLabels, Site, StructuralLabelsLess>{};
    for (const auto& chain : model.chains) {
        for (const auto& residue : chain.residues) {
            for (const auto& atom : residue.atoms) {
                const auto [site, added] =
                    by_labels.try_emplace(selected_labels({.chain = std::addressof(chain),
                                                           .residue = std::addressof(residue),
                                                           .atom = std::addressof(atom)}),
                                          Site{.atom = std::addressof(atom), .id = next_id});
                site->second.repeated = site->second.repeated || !added;
                ++next_id;
            }
        }
    }
    auto result = std::vector<std::string>{};
    result.reserve(conformer.sites.size());
    for (std::size_t index = 0; index < conformer.sites.size(); ++index) {
        const auto& labels = conformer.sites[index].structural_labels;
        if (labels == nullptr) {
            throw std::invalid_argument{"PDB annotation requires the labels of the PDB input"};
        }
        auto identity = *labels;
        identity.author.sequence = normalized_sequence(labels->author.sequence);
        const auto found = by_labels.find(identity);
        if (found == by_labels.end()) {
            throw std::invalid_argument{"PDB source does not match the calculation input"};
        }
        const auto& site = found->second;
        if (site.repeated) {
            throw std::invalid_argument{std::format(
                "PDB annotation cannot tell atom {} of residue {} {} {} from another atom with "
                "the same labels",
                site.atom->name, identity.author.residue.value_or(""),
                identity.author.chain.value_or(""), identity.author.sequence.value_or(""))};
        }
        const auto& position = molecule.conformer(conformer_index)[index];
        if (site.atom->element.atomic_number() != molecule.atom(index).atomic_number() ||
            site.atom->pos.x != position.x || site.atom->pos.y != position.y ||
            site.atom->pos.z != position.z) {
            throw std::invalid_argument{"PDB source does not match the calculation input"};
        }
        result.push_back(std::to_string(site.id));
    }
    return result;
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
    // The block is named as the molecule, so that reading the file gives the same record ID.
    structure.name = record.identity.record_id.display_string();
    if (structure.name.empty()) {
        structure.name = "structure";
    }

    auto result = ConvertedPdb{};
    auto next_id = std::size_t{1};
    for (std::size_t index = 0; index < metadata.conformers.size(); ++index) {
        result.atom_site_ids.push_back(calculated_site_ids(
            structure.models[index], next_id, metadata.conformers[index], record.molecule, index));
    }
    result.document = ::gemmi::make_mmcif_document(structure);
    for (auto& character : result.document.blocks.front().name) {
        if (std::isgraph(static_cast<unsigned char>(character)) == 0) {
            character = '_';
        }
    }
    return result;
}

} // namespace chargefw::io
