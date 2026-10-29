#include <chargefw/adapters/gemmi/pdb_input.h>

#include "bonds.h"
#include "selection.h"
#include "structure_import.h"

#include <chargefw/core/bond.h>

#include <gemmi/pdb.hpp>
#include <gemmi/polyheur.hpp>

#include <algorithm>
#include <cctype>
#include <charconv>
#include <cstdlib>
#include <istream>
#include <iterator>
#include <map>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

namespace chargefw::adapters::gemmi::pdb_input {
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

[[nodiscard]] auto parse_source(const std::string& contents) -> PdbSource {
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
        current->sites.push_back(
            PdbSourceSite{.reference = SourceAtomReference{.position = source_position++,
                                                           .id = field(view, 6, 5),
                                                           .structural_labels = std::move(labels)},
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

[[nodiscard]] auto source_models(const ::gemmi::Structure& structure, const InputOptions options,
                                 const std::span<const selection::SelectedModel> selected_models,
                                 const std::vector<PdbSourceModel>& source)
    -> std::vector<structure_import::SourceModelMapping> {
    const auto retained_count =
        options.conformers == ConformerSelection::all ? structure.models.size() : std::size_t{1};
    if (source.size() != structure.models.size()) {
        throw std::runtime_error{"PDB source mapping does not match parsed models"};
    }
    auto result = std::vector<structure_import::SourceModelMapping>{};
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
        result.push_back(structure_import::SourceModelMapping{
            .conformer = SourceConformerReference{
                .position = model_index, .id = source[model_index].id, .sites = std::move(sites)}});
    }
    return result;
}

} // namespace

PdbReader::PdbReader(std::istream& input, std::string source,
                     const ::chargefw::adapters::gemmi::InputOptions options)
    : options_{options} {
    std::string contents{std::istreambuf_iterator<char>{input}, std::istreambuf_iterator<char>{}};
    if (input.bad()) {
        throw std::runtime_error{"failed to read PDB input"};
    }

    auto structure = ::gemmi::read_pdb_string(contents, source);
    if (structure.models.empty()) {
        throw std::runtime_error{"structural input contains no models"};
    }
    ::gemmi::add_entity_types(structure, false);
    const auto name = structure.name.empty() ? source : structure.name;
    auto selected_models =
        selection::select_models(structure, options.selection, options.conformers);
    const auto parsed_source = parse_source(contents);
    auto mappings = source_models(structure, options, selected_models, parsed_source.models);
    for (std::size_t index = 0; index < selected_models.size(); ++index) {
        selected_models[index].normalize_component_names(mappings[index].conformer.sites);
    }
    auto explicit_bonds = bonds::ExplicitBonds{};
    if (options.bond_strategy == BondStrategy::explicit_bonds ||
        options.bond_strategy == BondStrategy::hybrid) {
        explicit_bonds = bonds::explicit_pdb(structure, selected_models.front());
    }
    record_ = structure_import::make_record(
        structure, selected_models,
        MoleculeRecordIdentity{.source = std::move(source), .record_index = 0, .record_id = name},
        options.selection, options.bond_strategy, options.conformers,
        std::move(explicit_bonds.bonds), name, MolecularSourceFormat::pdb, std::move(mappings),
        parsed_source.connectivity);
    if (explicit_bonds.skipped_symmetry_connections > 0) {
        record_->diagnostics.push_back(
            bonds::symmetry_connections_diagnostic(explicit_bonds.skipped_symmetry_connections));
    }
}

auto PdbReader::next() -> std::optional<ImportedMoleculeRecord> {
    return std::exchange(record_, std::nullopt);
}

auto PdbReader::options() const noexcept -> ::chargefw::adapters::gemmi::InputOptions {
    return options_;
}

} // namespace chargefw::adapters::gemmi::pdb_input
