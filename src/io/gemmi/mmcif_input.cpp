#include <chargefw/io/mmcif_input.h>
#include <chargefw/io/parse_error.h>

#include "bonds.h"
#include "io/imported_records.h"
#include "io/parse_context.h"
#include "io/text_encoding.h"
#include "mmcif_labels.h"
#include "selection.h"
#include "structure_import.h"

#include <chargefw/core/bond.h>

#include <gemmi/cif.hpp>
#include <gemmi/mmcif.hpp>
#include <gemmi/polyheur.hpp>

#include <algorithm>
#include <charconv>
#include <cstddef>
#include <istream>
#include <iterator>
#include <limits>
#include <memory>
#include <new>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <unordered_set>
#include <utility>
#include <vector>

namespace chargefw::io {
namespace {

// Replaces the validated _atom_site.id values with their one-based positions for Gemmi while it
// builds the structure, and restores the source IDs when destroyed. Renumbering in place avoids
// copying the whole block.
class ParserAtomSiteIds {
  public:
    explicit ParserAtomSiteIds(::gemmi::cif::Block& block)
        : atom_sites_{block.find("_atom_site.", {"id"})} {
        if (atom_sites_.length() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
            throw ParseError{
                ParseErrorKind::unsupported, "mmCIF contains too many atom sites for Gemmi", {}};
        }
        auto ids = std::unordered_set<std::string>{};
        ids.reserve(atom_sites_.length());
        for (auto row : atom_sites_) {
            if (::gemmi::cif::is_null(row[0])) {
                throw ParseError{ParseErrorKind::malformed,
                                 "mmCIF _atom_site.id must not be missing or unknown",
                                 {}};
            }
            auto source_id = ::gemmi::cif::as_string(row[0]);
            if (source_id.empty()) {
                throw ParseError{ParseErrorKind::malformed,
                                 "mmCIF _atom_site.id must not be missing or unknown",
                                 {}};
            }
            if (!ids.insert(source_id).second) {
                throw ParseError{ParseErrorKind::malformed,
                                 "mmCIF _atom_site.id must be unique; duplicate value '" +
                                     source_id + "'",
                                 {}};
            }
        }
        source_ids_.reserve(atom_sites_.length());
        try {
            auto position = std::size_t{0};
            for (auto row : atom_sites_) {
                auto parser_id = std::to_string(++position);
                source_ids_.push_back(std::exchange(row[0], std::move(parser_id)));
            }
        } catch (...) {
            // The destructor does not run when construction fails.
            restore();
            throw;
        }
    }

    ParserAtomSiteIds(const ParserAtomSiteIds&) = delete;
    ParserAtomSiteIds(ParserAtomSiteIds&&) = delete;
    auto operator=(const ParserAtomSiteIds&) -> ParserAtomSiteIds& = delete;
    auto operator=(ParserAtomSiteIds&&) -> ParserAtomSiteIds& = delete;

    ~ParserAtomSiteIds() {
        restore();
    }

  private:
    // Puts back the source IDs saved so far; later rows were not changed.
    auto restore() noexcept -> void {
        for (std::size_t index = 0; index < source_ids_.size(); ++index) {
            atom_sites_[static_cast<int>(index)][0] = std::move(source_ids_[index]);
        }
    }

    ::gemmi::cif::Table atom_sites_;
    std::vector<std::string> source_ids_;
};

[[nodiscard]] auto
make_source_mappings(::gemmi::cif::Block& block,
                     const std::span<const gemmi::selection::SelectedModel> selected_models)
    -> std::vector<gemmi::structure_import::SourceModelMapping> {
    // Gemmi built the structure while each _atom_site.id held its row position plus one, so an
    // atom's serial number names its row. Only the rows of selected atoms are decoded.
    auto atom_sites = gemmi::mmcif_labels::source_atom_sites(block);
    auto result = std::vector<gemmi::structure_import::SourceModelMapping>{};
    result.reserve(selected_models.size());
    for (std::size_t model_index = 0; model_index < selected_models.size(); ++model_index) {
        const auto& selected = selected_models[model_index];
        auto mapped = std::vector<SourceAtomReference>{};
        mapped.reserve(selected.atoms().size());
        auto model_id = std::optional<std::string>{};
        auto model_id_initialized = false;
        for (const auto* atom : selected.atoms()) {
            if (atom->serial <= 0 || static_cast<std::size_t>(atom->serial) > atom_sites.length()) {
                throw std::logic_error{"mmCIF source mapping does not match selected atoms"};
            }
            const auto row = atom_sites[atom->serial - 1];
            if (!model_id_initialized) {
                model_id =
                    gemmi::mmcif_labels::source_value(row, gemmi::mmcif_labels::model_id_column);
                model_id_initialized = true;
            }
            mapped.push_back(SourceAtomReference{
                .position = static_cast<std::size_t>(atom->serial) - 1,
                .id = ::gemmi::cif::as_string(row[0]),
                .structural_labels = std::make_shared<const SourceStructuralLabels>(
                    gemmi::mmcif_labels::decode_source_labels(row))});
        }
        result.push_back(gemmi::structure_import::SourceModelMapping{
            .conformer = SourceConformerReference{
                .position = model_index, .id = std::move(model_id), .sites = std::move(mapped)}});
    }
    return result;
}

[[nodiscard]] auto source_connectivity(const ::gemmi::cif::Block& block) -> SourceConnectivity {
    return block.has_mmcif_category("_chem_comp_bond.") || block.has_mmcif_category("_struct_conn.")
               ? SourceConnectivity::present
               : SourceConnectivity::absent;
}

} // namespace

struct MmcifReader::Impl {
    Impl(std::istream& input, std::string source_name, ::chargefw::io::StructureOptions options);
    auto next() -> std::optional<ImportedMoleculeRecord>;
    auto read_block(::gemmi::cif::Block& block, std::size_t record_index) -> ImportedMoleculeRecord;

    ::gemmi::cif::Document document_;
    // Kept by the records for annotated output: a tenth of the memory of the parsed document.
    std::shared_ptr<const std::string> text_;
    std::string source_;
    ::chargefw::io::StructureOptions options_;
    std::size_t block_index_ = 0;
    std::size_t record_index_ = 0;
};

MmcifReader::MmcifReader(std::istream& input, std::string source_name,
                         const ::chargefw::io::StructureOptions options)
    : impl_{std::make_unique<Impl>(input, std::move(source_name), options)} {}

MmcifReader::MmcifReader(MmcifReader&&) noexcept = default;
auto MmcifReader::operator=(MmcifReader&&) noexcept -> MmcifReader& = default;
MmcifReader::~MmcifReader() = default;

auto MmcifReader::next() -> std::optional<ImportedMoleculeRecord> {
    return impl_ ? impl_->next() : std::nullopt;
}

namespace {

// A formal charge Gemmi cannot hold would be read as another charge, so it is refused first.
auto require_readable_formal_charges(::gemmi::cif::Block& block) -> void {
    for (const auto row : block.find("_atom_site.", {"pdbx_formal_charge"})) {
        if (!gemmi::mmcif_labels::formal_charge(row[0]).has_value()) {
            throw ParseError{ParseErrorKind::malformed,
                             "mmCIF _atom_site.pdbx_formal_charge '" +
                                 ::gemmi::cif::as_string(row[0]) + "' is not a formal charge",
                             {}};
        }
    }
}

// Gemmi reports an error in the text as "source:line:column(byte): reason", as
// "source:line in data_block: reason", or without a line as "source: reason"; the line and the
// reason are taken from it, or the whole message is the reason.
[[nodiscard]] auto syntax_error(std::string_view message, const std::string& source) -> ParseError {
    auto location = ParseLocation{.source_name = source};
    const auto prefix = source + ":";
    if (message.starts_with(prefix)) {
        const auto rest = message.substr(prefix.size());
        auto line = std::size_t{};
        const auto [end, error] = std::from_chars(rest.data(), rest.data() + rest.size(), line);
        const auto reason = rest.find(": ", static_cast<std::size_t>(end - rest.data()));
        if (error == std::errc{} && reason != std::string_view::npos) {
            location.line = line;
            message = rest.substr(reason + 2);
        } else if (rest.starts_with(' ')) {
            message = rest.substr(1);
        }
    }
    return ParseError{ParseErrorKind::malformed, std::string{message}, std::move(location)};
}

} // namespace

MmcifReader::Impl::Impl(std::istream& input, std::string source_name,
                        const ::chargefw::io::StructureOptions options)
    : source_{std::move(source_name)}, options_{options} {
    // Parse from memory: a buffered stream parser rejects tokens longer than its buffer.
    io::detail::require_readable(input);
    auto text =
        std::string{std::istreambuf_iterator<char>{input}, std::istreambuf_iterator<char>{}};
    io::detail::require_readable(input);
    if (const auto line = io::detail::first_non_ascii_line(text)) {
        throw ParseError{ParseErrorKind::invalid_encoding,
                         "non-ASCII character; mmCIF files must be ASCII text",
                         {.source_name = source_, .line = line}};
    }
    try {
        document_ = ::gemmi::cif::read_memory(text.data(), text.size(), source_.c_str());
    } catch (const std::bad_alloc&) {
        throw;
    } catch (const std::exception& error) {
        throw syntax_error(error.what(), source_);
    }
    text_ = std::make_shared<const std::string>(std::move(text));

    const auto has_coordinates = [](const ::gemmi::cif::Block& block) -> bool {
        return block.has_tag("_atom_site.id");
    };
    if (std::ranges::none_of(document_.blocks, has_coordinates)) {
        throw ParseError{ParseErrorKind::no_records,
                         "mmCIF input contains no coordinate data blocks",
                         {.source_name = source_}};
    }
}

auto MmcifReader::Impl::next() -> std::optional<ImportedMoleculeRecord> {
    while (block_index_ < document_.blocks.size()) {
        auto& block = document_.blocks[block_index_++];
        if (!block.has_tag("_atom_site.id")) {
            continue;
        }
        const auto record_index = record_index_++;
        try {
            return read_block(block, record_index);
        } catch (...) {
            io::detail::rethrow_located({.source_name = source_,
                                         .line = std::nullopt,
                                         .record_index = record_index,
                                         .record_id = block.name});
        }
    }

    return std::nullopt;
}

auto MmcifReader::Impl::read_block(::gemmi::cif::Block& block, const std::size_t record_index)
    -> ImportedMoleculeRecord {
    auto structure = io::detail::parsed([&block] {
        require_readable_formal_charges(block);
        const auto parser_ids = ParserAtomSiteIds{block};
        return ::gemmi::make_structure_from_block(block);
    });
    if (structure.models.empty()) {
        throw ParseError{ParseErrorKind::malformed, "mmCIF block has no atom sites", {}};
    }
    ::gemmi::add_entity_types(structure, false);
    // Atoms are told apart by their label name where the file has one, as residue templates match
    // them; Gemmi names an atom by its author name. An atom's serial number names its row, as in
    // make_source_mappings().
    auto atom_sites = gemmi::mmcif_labels::source_atom_sites(block);
    const auto atom_name = [&atom_sites](const ::gemmi::Atom& atom) {
        const auto label = gemmi::mmcif_labels::source_value(atom_sites[atom.serial - 1], 1);
        return std::string{io::detail::known_label(label).value_or(atom.name)};
    };
    auto selected_models = gemmi::selection::select_models(structure, options_.selection,
                                                           options_.conformers, atom_name);
    auto source_models = make_source_mappings(block, selected_models);
    for (std::size_t index = 0; index < selected_models.size(); ++index) {
        selected_models[index].normalize_component_names(source_models[index].conformer.sites);
    }
    auto explicit_bonds = gemmi::bonds::ExplicitBonds{};
    if (options_.bonds == BondStrategy::explicit_bonds || options_.bonds == BondStrategy::hybrid) {
        explicit_bonds =
            gemmi::bonds::explicit_mmcif(structure, block, selected_models.front(), options_.bonds);
    }
    auto record = gemmi::structure_import::make_record(
        structure, selected_models,
        MoleculeRecordIdentity{
            .source_name = source_, .record_index = record_index, .record_id = block.name},
        options_.selection, options_.bonds, options_.conformers, std::move(explicit_bonds.bonds),
        structure.name.empty() ? block.name : structure.name, MolecularSourceFormat::mmcif,
        std::move(source_models), source_connectivity(block));
    record.import_metadata->source_text = text_;
    if (!explicit_bonds.aromatic_bonds.empty()) {
        record.import_warnings.push_back(
            io::aromatic_bonds_issue(std::move(explicit_bonds.aromatic_bonds)));
    }
    if (explicit_bonds.skipped_symmetry_connections > 0) {
        record.import_warnings.push_back(
            gemmi::bonds::symmetry_connections_issue(explicit_bonds.skipped_symmetry_connections));
    }
    // Either category, as the annotating writer refuses a block that has either.
    if (block.has_mmcif_category("_sb_ncbr_partial_atomic_charges.") ||
        block.has_mmcif_category("_sb_ncbr_partial_atomic_charges_meta.")) {
        record.import_warnings.push_back(ImportIssue{
            .kind = ImportIssueKind::partial_charges_ignored,
            .message = "mmCIF partial charges were ignored; they are not formal charges and "
                       "were not used for calculation",
            .line = std::nullopt});
    }
    return record;
}

} // namespace chargefw::io
