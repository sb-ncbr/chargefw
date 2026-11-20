#include "pdb_conversion.h"

#include "pdb_labels.h"

#include <chargefw/core/issue_location.h>
#include <chargefw/core/molecule.h>
#include <chargefw/io/molecule_record.h>
#include <chargefw/io/output_error.h>

#include <gemmi/align.hpp>
#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include <gemmi/polyheur.hpp>
#include <gemmi/to_mmcif.hpp>

#include <cctype>
#include <cstddef>
#include <format>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace chargefw::io {
namespace {

using gemmi::pdb_labels::normalized_sequence;
using gemmi::pdb_labels::selected_labels;
using gemmi::pdb_labels::StructuralLabelsLess;

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
            throw OutputError{
                OutputErrorKind::unrepresentable,
                std::format("PDB annotation cannot tell atom {} of residue {} {} {} from another "
                            "atom with the same labels",
                            site.atom->name, identity.author.residue.value_or(""),
                            identity.author.chain.value_or(""),
                            identity.author.sequence.value_or("")),
                0,
                core::IssueLocation{.conformer_index = conformer_index,
                                    .atom_index = index,
                                    .atom_indices = {index}}};
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
