#include "selection.h"

#include <algorithm>
#include <chargefw/io/parse_error.h>
#include <iterator>
#include <memory>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>

namespace chargefw::io::gemmi::selection {
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

SelectedModel::SelectedModel(const ::gemmi::Model& model, const RecordSelection selection,
                             const AtomName& atom_name) {
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
            // Of alternative residues at one position, such as SER and PRO at 22 with alternate
            // locations A and C, read the first. Gemmi groups adjacent residues by sequence number
            // and insertion code alone, and files also repeat a number for residues that are all
            // there. A residue is an alternative of the one before it only in the same PDB segment
            // and mmCIF label chain, and when all its atoms have an alternate location. Residues
            // without a number are not alternatives.
            const auto alternative = [&] {
                if (!residue.seqid.num.has_value() || chain.is_first_in_group(residue)) {
                    return false;
                }
                const auto& previous = *std::prev(std::addressof(residue));
                return previous.segment == residue.segment &&
                       previous.subchain == residue.subchain &&
                       std::ranges::all_of(residue.atoms,
                                           [](const auto& atom) { return atom.has_altloc(); });
            }();
            if (alternative || !include_residue(residue, selection)) {
                continue;
            }

            SelectedResidue selected{.residue = std::addressof(residue),
                                     .chain_name = chain.name,
                                     .component_name = residue.name,
                                     .atoms = {}};
            selected.atoms.reserve(residue.atoms.size());
            // An atom name repeats in a residue only for alternate locations, of which the first
            // in source order is read. Copies that no alternate location tells apart are not
            // alternatives: the file repeats an atom or gave several residues the same labels.
            auto seen_altlocs = std::unordered_map<std::string, std::string>{};
            seen_altlocs.reserve(residue.atoms.size());
            for (std::size_t index = 0; index < residue.atoms.size(); ++index) {
                const auto& atom = residue.atoms[index];
                auto name = atom_name ? atom_name(atom) : atom.name;
                auto& altlocs = seen_altlocs[name];
                if (altlocs.contains(atom.altloc)) {
                    const auto icode = residue.seqid.has_icode() ? std::string{residue.seqid.icode}
                                                                 : std::string{};
                    throw ParseError{ParseErrorKind::malformed,
                                     "residue " +
                                         detail::residue_description(residue.name, chain.name,
                                                                     residue.seqid.num.str(),
                                                                     icode) +
                                         " has more than one atom named " + name +
                                         " and no alternate location to tell them apart; the atom "
                                         "is repeated, or "
                                         "several residues have this name, chain, and number",
                                     {}};
                }
                const auto alternate = !altlocs.empty();
                altlocs.push_back(atom.altloc);
                if (alternate) {
                    continue;
                }

                atoms_.push_back(std::addressof(atom));
                sites_.push_back(::gemmi::const_CRA{.chain = std::addressof(chain),
                                                    .residue = std::addressof(residue),
                                                    .atom = std::addressof(atom)});
                selected.atoms.push_back({std::move(name), atom_index});
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
        throw std::logic_error{"structural source mapping does not match selected atoms"};
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
                throw ParseError{ParseErrorKind::malformed,
                                 "selected structural component has an empty canonical name",
                                 {}};
            }
            if (component_name.empty()) {
                component_name = canonical_component;
            } else if (component_name != canonical_component) {
                throw ParseError{ParseErrorKind::malformed,
                                 "selected structural residue has inconsistent component names",
                                 {}};
            }
            if (!names.emplace(canonical_atom).second) {
                throw ParseError{ParseErrorKind::malformed,
                                 "selected structural component has ambiguous canonical atom names",
                                 {}};
            }
            atom.name = canonical_atom;
        }
        residue.component_name = std::move(component_name);
    }
}

auto select_models(const ::gemmi::Structure& structure, const RecordSelection selection,
                   const ConformerSelection conformers, const AtomName& atom_name)
    -> std::vector<SelectedModel> {
    const auto count =
        conformers == ConformerSelection::all ? structure.models.size() : std::size_t{1};
    auto result = std::vector<SelectedModel>{};
    result.reserve(count);
    for (std::size_t index = 0; index < count; ++index) {
        result.emplace_back(structure.models[index], selection, atom_name);
    }
    return result;
}

} // namespace chargefw::io::gemmi::selection
