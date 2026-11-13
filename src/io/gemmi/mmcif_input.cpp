#include <chargefw/io/mmcif_input.h>

#include "bonds.h"
#include "io/imported_records.h"
#include "io/text_encoding.h"
#include "mmcif_labels.h"
#include "selection.h"
#include "structure_import.h"

#include <chargefw/core/bond.h>

#include <gemmi/cif.hpp>
#include <gemmi/mmcif.hpp>
#include <gemmi/polyheur.hpp>

#include <algorithm>
#include <cstddef>
#include <istream>
#include <iterator>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
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
            throw std::runtime_error{"mmCIF contains too many atom sites for Gemmi"};
        }
        auto ids = std::unordered_set<std::string>{};
        ids.reserve(atom_sites_.length());
        for (auto row : atom_sites_) {
            if (::gemmi::cif::is_null(row[0])) {
                throw std::runtime_error{"mmCIF _atom_site.id must not be missing or unknown"};
            }
            auto source_id = ::gemmi::cif::as_string(row[0]);
            if (source_id.empty()) {
                throw std::runtime_error{"mmCIF _atom_site.id must not be missing or unknown"};
            }
            if (!ids.insert(source_id).second) {
                throw std::runtime_error{"mmCIF _atom_site.id must be unique; duplicate value '" +
                                         source_id + "'"};
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
                throw std::runtime_error{"mmCIF source mapping does not match selected atoms"};
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

    ::gemmi::cif::Document document_;
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

MmcifReader::Impl::Impl(std::istream& input, std::string source_name,
                        const ::chargefw::io::StructureOptions options)
    : source_{std::move(source_name)}, options_{options} {
    // Parse from memory: a buffered stream parser rejects tokens longer than its buffer.
    const std::string text{std::istreambuf_iterator<char>{input}, std::istreambuf_iterator<char>{}};
    if (input.bad()) {
        throw std::runtime_error{"failed to read mmCIF input"};
    }
    if (const auto line = io::detail::first_non_ascii_line(text)) {
        throw std::runtime_error{"line " + std::to_string(*line) +
                                 " contains a non-ASCII character; mmCIF files must be ASCII text"};
    }
    document_ = ::gemmi::cif::read_memory(text.data(), text.size(), source_.c_str());

    const auto has_coordinates = [](const ::gemmi::cif::Block& block) -> bool {
        return block.has_tag("_atom_site.id");
    };
    if (std::ranges::none_of(document_.blocks, has_coordinates)) {
        throw std::runtime_error{"mmCIF input contains no coordinate data blocks"};
    }
}

auto MmcifReader::Impl::next() -> std::optional<ImportedMoleculeRecord> {
    while (block_index_ < document_.blocks.size()) {
        auto& block = document_.blocks[block_index_++];
        if (!block.has_tag("_atom_site.id")) {
            continue;
        }

        auto structure = [&block] {
            const auto parser_ids = ParserAtomSiteIds{block};
            return ::gemmi::make_structure_from_block(block);
        }();
        if (structure.models.empty()) {
            throw std::runtime_error{"structural input contains no models"};
        }
        ::gemmi::add_entity_types(structure, false);
        auto selected_models =
            gemmi::selection::select_models(structure, options_.selection, options_.conformers);
        auto source_models = make_source_mappings(block, selected_models);
        for (std::size_t index = 0; index < selected_models.size(); ++index) {
            selected_models[index].normalize_component_names(source_models[index].conformer.sites);
        }
        auto explicit_bonds = gemmi::bonds::ExplicitBonds{};
        if (options_.bonds == BondStrategy::explicit_bonds ||
            options_.bonds == BondStrategy::hybrid) {
            explicit_bonds =
                gemmi::bonds::explicit_mmcif(structure, block, selected_models.front());
        }
        const auto current_record_index = record_index_++;
        auto record = gemmi::structure_import::make_record(
            structure, selected_models,
            MoleculeRecordIdentity{.source_name = source_,
                                   .record_index = current_record_index,
                                   .record_id = block.name},
            options_.selection, options_.bonds, options_.conformers,
            std::move(explicit_bonds.bonds), structure.name.empty() ? block.name : structure.name,
            MolecularSourceFormat::mmcif, std::move(source_models), source_connectivity(block));
        if (!explicit_bonds.aromatic_bonds.empty()) {
            record.import_warnings.push_back(
                io::aromatic_bonds_issue(std::move(explicit_bonds.aromatic_bonds)));
        }
        if (explicit_bonds.skipped_symmetry_connections > 0) {
            record.import_warnings.push_back(gemmi::bonds::symmetry_connections_issue(
                explicit_bonds.skipped_symmetry_connections));
        }
        if (block.has_mmcif_category("_sb_ncbr_partial_atomic_charges.")) {
            record.import_warnings.push_back(ImportIssue{
                .kind = ImportIssueKind::partial_charges_ignored,
                .message = "mmCIF partial charges were ignored; they are not formal charges and "
                           "were not used for calculation.",
                .line = std::nullopt});
        }
        return record;
    }

    return std::nullopt;
}

} // namespace chargefw::io
