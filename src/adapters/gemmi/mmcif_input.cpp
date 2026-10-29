#include <chargefw/adapters/gemmi/mmcif_input.h>

#include "bonds.h"
#include "mmcif_labels.h"
#include "selection.h"
#include "structure_import.h"

#include <chargefw/core/bond.h>

#include <gemmi/cif.hpp>
#include <gemmi/mmcif.hpp>
#include <gemmi/polyheur.hpp>

#include <algorithm>
#include <istream>
#include <iterator>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>

namespace chargefw::adapters::gemmi::mmcif_input {
namespace {

[[nodiscard]] auto parser_block(const ::gemmi::cif::Block& source) -> ::gemmi::cif::Block {
    auto result = source;
    auto atom_sites = result.find("_atom_site.", {"id"});
    if (atom_sites.length() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::runtime_error{"mmCIF contains too many atom sites for Gemmi"};
    }
    auto ids = std::unordered_set<std::string>{};
    ids.reserve(atom_sites.length());
    auto position = std::size_t{0};
    for (auto row : atom_sites) {
        if (::gemmi::cif::is_null(row[0])) {
            throw std::runtime_error{"mmCIF _atom_site.id must not be missing or unknown"};
        }
        const auto source_id = ::gemmi::cif::as_string(row[0]);
        if (source_id.empty()) {
            throw std::runtime_error{"mmCIF _atom_site.id must not be missing or unknown"};
        }
        if (!ids.insert(source_id).second) {
            throw std::runtime_error{"mmCIF _atom_site.id must be unique; duplicate value '" +
                                     source_id + "'"};
        }
        row[0] = std::to_string(++position);
    }
    return result;
}

struct MmcifSourceSite {
    std::optional<std::string> model_id;
    SourceAtomReference reference;
};

[[nodiscard]] auto source_sites(::gemmi::cif::Block& block) -> std::vector<MmcifSourceSite> {
    auto atom_sites = mmcif_labels::source_atom_sites(block);
    auto result = std::vector<MmcifSourceSite>{};
    result.reserve(atom_sites.length());
    auto position = std::size_t{0};
    for (const auto row : atom_sites) {
        const auto id = ::gemmi::cif::as_string(row[0]);
        const auto model_id = mmcif_labels::source_value(row, mmcif_labels::model_id_column);
        result.push_back(MmcifSourceSite{
            .model_id = model_id,
            .reference =
                SourceAtomReference{.position = position++,
                                    .id = id,
                                    .structural_labels = mmcif_labels::decode_source_labels(row)}});
    }
    return result;
}

[[nodiscard]] auto
make_source_mappings(::gemmi::cif::Block& block,
                     const std::span<const selection::SelectedModel> selected_models)
    -> std::vector<structure_import::SourceModelMapping> {
    const auto sites = source_sites(block);
    auto result = std::vector<structure_import::SourceModelMapping>{};
    result.reserve(selected_models.size());
    for (std::size_t model_index = 0; model_index < selected_models.size(); ++model_index) {
        const auto& selected = selected_models[model_index];
        auto mapped = std::vector<SourceAtomReference>{};
        mapped.reserve(selected.atoms().size());
        auto model_id = std::optional<std::string>{};
        auto model_id_initialized = false;
        for (const auto* atom : selected.atoms()) {
            if (atom->serial <= 0 || static_cast<std::size_t>(atom->serial) > sites.size()) {
                throw std::runtime_error{"mmCIF source mapping does not match selected atoms"};
            }
            const auto& source = sites[static_cast<std::size_t>(atom->serial) - 1];
            if (!model_id_initialized) {
                model_id = source.model_id;
                model_id_initialized = true;
            }
            mapped.push_back(source.reference);
        }
        result.push_back(structure_import::SourceModelMapping{
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

MmcifReader::MmcifReader(std::istream& input, std::string source,
                         const ::chargefw::adapters::gemmi::InputOptions options)
    : source_{std::move(source)}, options_{options} {
    // Parse from memory: a buffered stream parser rejects tokens longer than its buffer.
    const std::string text{std::istreambuf_iterator<char>{input}, std::istreambuf_iterator<char>{}};
    if (input.bad()) {
        throw std::runtime_error{"failed to read mmCIF input"};
    }
    document_ = ::gemmi::cif::read_memory(text.data(), text.size(), source_.c_str());

    const auto has_coordinates = [](const ::gemmi::cif::Block& block) -> bool {
        return block.has_tag("_atom_site.id");
    };
    if (std::ranges::none_of(document_.blocks, has_coordinates)) {
        throw std::runtime_error{"mmCIF input contains no coordinate data blocks"};
    }
}

auto MmcifReader::next() -> std::optional<ImportedMoleculeRecord> {
    while (block_index_ < document_.blocks.size()) {
        auto& block = document_.blocks[block_index_++];
        if (!block.has_tag("_atom_site.id")) {
            continue;
        }

        auto normalized = parser_block(block);
        auto structure = ::gemmi::make_structure_from_block(normalized);
        if (structure.models.empty()) {
            throw std::runtime_error{"structural input contains no models"};
        }
        ::gemmi::add_entity_types(structure, false);
        auto selected_models =
            selection::select_models(structure, options_.selection, options_.conformers);
        auto source_models = make_source_mappings(block, selected_models);
        for (std::size_t index = 0; index < selected_models.size(); ++index) {
            selected_models[index].normalize_component_names(source_models[index].conformer.sites);
        }
        auto explicit_bonds = std::vector<core::Bond>{};
        if (options_.bond_strategy == BondStrategy::explicit_bonds ||
            options_.bond_strategy == BondStrategy::hybrid) {
            explicit_bonds = bonds::explicit_mmcif(structure, block, selected_models.front());
        }
        const auto current_record_index = record_index_++;
        auto record = structure_import::make_record(
            structure, selected_models,
            MoleculeRecordIdentity{
                .source = source_, .record_index = current_record_index, .record_id = block.name},
            options_.selection, options_.bond_strategy, options_.conformers,
            std::move(explicit_bonds), structure.name.empty() ? block.name : structure.name,
            MolecularSourceFormat::mmcif, std::move(source_models), source_connectivity(block));
        if (block.has_mmcif_category("_sb_ncbr_partial_atomic_charges.")) {
            record.diagnostics.push_back(MoleculeRecordDiagnostic{
                .code = "partial_charges_ignored",
                .message = "mmCIF partial charges were ignored; they are not formal charges and "
                           "were not used for calculation.",
                .line = std::nullopt});
        }
        return record;
    }

    return std::nullopt;
}

auto MmcifReader::options() const noexcept -> ::chargefw::adapters::gemmi::InputOptions {
    return options_;
}

} // namespace chargefw::adapters::gemmi::mmcif_input
