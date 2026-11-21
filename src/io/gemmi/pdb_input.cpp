#include <chargefw/io/parse_error.h>
#include <chargefw/io/pdb_input.h>

#include "bonds.h"
#include "io/parse_context.h"
#include "io/text_encoding.h"
#include "pdb_labels.h"
#include "selection.h"
#include "structure_import.h"

#include <chargefw/core/bond.h>

#include <gemmi/calculate.hpp>
#include <gemmi/pdb.hpp>
#include <gemmi/polyheur.hpp>

#include <algorithm>
#include <cstddef>
#include <istream>
#include <iterator>
#include <map>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace chargefw::io {
namespace {

using gemmi::pdb_labels::normalized_sequence;
using gemmi::pdb_labels::optional_text;
using gemmi::pdb_labels::selected_labels;
using gemmi::pdb_labels::StructuralLabelsLess;

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

[[nodiscard]] auto
source_models(const ::gemmi::Structure& structure, const StructureOptions& options,
              const std::span<const gemmi::selection::SelectedModel> selected_models,
              const std::vector<PdbSourceModel>& source)
    -> std::vector<gemmi::structure_import::SourceModelMapping> {
    const auto retained_count =
        options.conformers == ConformerSelection::all ? structure.models.size() : std::size_t{1};
    if (source.size() != structure.models.size()) {
        // Reachable from a file whose records Gemmi and ChargeFW's own pass read differently, such
        // as a residue number that is not a number, so it is a refusal and not an internal error.
        throw ParseError{
            ParseErrorKind::malformed, "PDB source mapping does not match parsed models", {}};
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
                throw ParseError{ParseErrorKind::malformed,
                                 "PDB source mapping does not match selected atoms",
                                 {}};
            }
            if (found->second.size() != 1) {
                throw ParseError{ParseErrorKind::malformed,
                                 "PDB source contains ambiguous selected atom identity",
                                 {}};
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

[[nodiscard]] auto read_record(std::string contents, std::string source_name,
                               const ::chargefw::io::StructureOptions& options)
    -> ImportedMoleculeRecord {
    if (const auto line = io::detail::first_non_ascii_line(contents)) {
        throw ParseError{ParseErrorKind::invalid_encoding,
                         "non-ASCII character; PDB files must be ASCII text",
                         {.source_name = source_name, .line = line}};
    }

    auto structure =
        io::detail::parsed([&] { return ::gemmi::read_pdb_string(contents, source_name); });
    // Gemmi gives a file without MODEL records one model, even when it has no atoms.
    if (::gemmi::count_atom_sites(structure) == 0) {
        throw ParseError{ParseErrorKind::no_records, "no PDB ATOM or HETATM record found", {}};
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
    auto record = gemmi::structure_import::make_record(
        structure, selected_models,
        MoleculeRecordIdentity{
            .source_name = std::move(source_name), .record_index = 0, .record_id = name},
        options.selection, options.bonds, options.components.get(), options.conformers,
        std::move(explicit_bonds.bonds), explicit_bonds.residues_with_file_bonds, name,
        MolecularSourceFormat::pdb, std::move(mappings), parsed_source.connectivity,
        std::move(explicit_bonds.unknown_order_indices));
    record.import_metadata->source_text = std::make_shared<const std::string>(std::move(contents));
    if (explicit_bonds.skipped_symmetry_connections > 0) {
        record.import_warnings.push_back(
            gemmi::bonds::symmetry_connections_issue(explicit_bonds.skipped_symmetry_connections));
    }
    return record;
}

} // namespace

PdbReader::PdbReader(std::istream& input, std::string source_name,
                     const ::chargefw::io::StructureOptions& options) {
    io::detail::require_readable(input);
    std::string contents{std::istreambuf_iterator<char>{input}, std::istreambuf_iterator<char>{}};
    io::detail::require_readable(input);
    const auto location = ParseLocation{.source_name = source_name};
    try {
        record_ = read_record(std::move(contents), std::move(source_name), options);
    } catch (...) {
        io::detail::rethrow_located(location);
    }
}

auto PdbReader::next() -> std::optional<ImportedMoleculeRecord> {
    return std::exchange(record_, std::nullopt);
}

} // namespace chargefw::io
