#include "selection.h"

#include <string_view>
#include <unordered_set>

namespace chargefw::adapters::gemmi::selection {
namespace {

auto include_residue(const ::gemmi::Residue& residue, const RecordSelection selection) -> bool {
    switch (selection) {
    case RecordSelection::all:
        return true;
    case RecordSelection::polymers_and_ligands:
        return !residue.is_water();
    case RecordSelection::polymers:
        return residue.entity_type == ::gemmi::EntityType::Polymer;
    }

    return false;
}

} // namespace

SelectedModel::SelectedModel(const ::gemmi::Model& model, const RecordSelection selection) {
    std::size_t atom_count = 0;
    for (const auto& chain : model.chains) {
        for (const auto& residue : chain.residues) {
            atom_count += residue.atoms.size();
        }
    }
    atoms_.reserve(atom_count);
    sites_.reserve(atom_count);
    residues_.reserve(model.chains.size());
    atom_indices_.reserve(atom_count);
    serial_indices_.reserve(atom_count);

    std::size_t atom_index = 0;
    for (const auto& chain : model.chains) {
        for (const auto& residue : chain.residues) {
            if (!include_residue(residue, selection)) {
                continue;
            }

            SelectedResidue selected{.residue = std::addressof(residue),
                                     .chain_name = chain.name,
                                     .component_name = residue.name,
                                     .atoms = {}};
            selected.atoms.reserve(residue.atoms.size());
            // Select one alternate location per atom name: the first source-order occurrence.
            auto seen_names = std::unordered_set<std::string_view>{};
            seen_names.reserve(residue.atoms.size());
            for (std::size_t index = 0; index < residue.atoms.size(); ++index) {
                const auto& atom = residue.atoms[index];
                if (!seen_names.emplace(atom.name).second) {
                    continue;
                }

                atoms_.push_back(std::addressof(atom));
                sites_.push_back(::gemmi::const_CRA{.chain = std::addressof(chain),
                                                    .residue = std::addressof(residue),
                                                    .atom = std::addressof(atom)});
                selected.atoms.push_back({atom.name, atom_index});
                atom_indices_.emplace(std::addressof(atom), atom_index);
                serial_indices_.emplace(atom.serial, atom_index++);
            }
            residues_.push_back(std::move(selected));
        }
    }
}

auto SelectedModel::atoms() const -> const std::vector<const ::gemmi::Atom*>& {
    return atoms_;
}

auto SelectedModel::sites() const -> const std::vector<::gemmi::const_CRA>& {
    return sites_;
}

auto SelectedModel::residues() const -> const std::vector<SelectedResidue>& {
    return residues_;
}

auto SelectedModel::atom_index(const ::gemmi::Atom* atom) const -> std::optional<std::size_t> {
    const auto found = atom_indices_.find(atom);
    if (found == atom_indices_.end()) {
        return std::nullopt;
    }
    return found->second;
}

auto SelectedModel::atom_index_by_serial(const int serial) const -> std::optional<std::size_t> {
    const auto found = serial_indices_.find(serial);
    if (found == serial_indices_.end()) {
        return std::nullopt;
    }
    return found->second;
}

void SelectedModel::normalize_component_names(
    const std::span<const SourceAtomReference> source_atoms) {
    if (source_atoms.size() != atoms_.size()) {
        throw std::runtime_error{"structural source mapping does not match selected atoms"};
    }
    for (auto& residue : residues_) {
        auto component_name = std::string{};
        auto names = std::unordered_set<std::string_view>{};
        names.reserve(residue.atoms.size());
        for (auto& atom : residue.atoms) {
            const auto* labels = source_atoms[atom.index].structural_labels.get();
            const auto canonical_component =
                detail::canonical_component_name(labels, residue.residue->name);
            const auto canonical_atom =
                detail::canonical_atom_name(labels, atoms_[atom.index]->name);
            if (canonical_component.empty() || canonical_atom.empty()) {
                throw std::runtime_error{
                    "selected structural component has an empty canonical name"};
            }
            if (component_name.empty()) {
                component_name = canonical_component;
            } else if (component_name != canonical_component) {
                throw std::runtime_error{
                    "selected structural residue has inconsistent component names"};
            }
            if (!names.emplace(canonical_atom).second) {
                throw std::runtime_error{
                    "selected structural component has ambiguous canonical atom names"};
            }
            atom.name = canonical_atom;
        }
        residue.component_name = std::move(component_name);
    }
}

auto select_models(const ::gemmi::Structure& structure, const RecordSelection selection,
                   const ConformerSelection conformers) -> std::vector<SelectedModel> {
    const auto count =
        conformers == ConformerSelection::all ? structure.models.size() : std::size_t{1};
    auto result = std::vector<SelectedModel>{};
    result.reserve(count);
    for (std::size_t index = 0; index < count; ++index) {
        result.emplace_back(structure.models[index], selection);
    }
    return result;
}

} // namespace chargefw::adapters::gemmi::selection
