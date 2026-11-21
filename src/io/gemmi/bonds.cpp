#include "bonds.h"

#include "component_templates.h"

#include "selection.h"
#include <chargefw/io/component_dictionary.h>

#include <chargefw/io/parse_error.h>
#include <gemmi/chemcomp.hpp>
#include <gemmi/cif.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdlib>
#include <iterator>
#include <map>
#include <memory>
#include <optional>
#include <ranges>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace chargefw::io::gemmi::bonds {
namespace {

// A bond order read from a file. Aromatic bonds are read as single bonds.
struct SourceOrder {
    core::BondOrder order;
    bool aromatic = false;
};

class BondAccumulator {
  public:
    auto add(const std::size_t first, const std::size_t second, const core::BondOrder order)
        -> std::optional<std::size_t> {
        return add(first, second, SourceOrder{.order = order});
    }

    // The first order given for a pair is kept. Return an index only for a newly inserted bond.
    auto add(const std::size_t first, const std::size_t second, const SourceOrder order)
        -> std::optional<std::size_t> {
        if (first == second) {
            return std::nullopt;
        }

        const auto key =
            BondKey{.first = std::min(first, second), .second = std::max(first, second)};
        if (seen_.insert(key).second) {
            bonds_.emplace_back(first, second, order.order);
            if (order.aromatic) {
                aromatic_bonds_.push_back(bonds_.size() - 1);
            }
            return bonds_.size() - 1;
        }
        return std::nullopt;
    }

    [[nodiscard]] auto aromatic_bonds() const -> const std::vector<std::size_t>& {
        return aromatic_bonds_;
    }

    [[nodiscard]] auto take() -> std::vector<core::Bond> {
        return std::move(bonds_);
    }

  private:
    struct BondKey {
        std::size_t first;
        std::size_t second;

        [[nodiscard]] auto operator==(const BondKey&) const -> bool = default;
    };

    struct BondKeyHash {
        [[nodiscard]] auto operator()(const BondKey& key) const noexcept -> std::size_t {
            return std::hash<std::size_t>{}(key.first) ^
                   (std::hash<std::size_t>{}(key.second) << 1U);
        }
    };

    std::vector<core::Bond> bonds_;
    std::unordered_set<BondKey, BondKeyHash> seen_;
    std::vector<std::size_t> aromatic_bonds_; // Indices in bonds_.
};

// The selected residues of a model by chain name and sequence ID, in model order, so that each
// connection partner is found without scanning every chain and residue.
class ResidueIndex {
  public:
    explicit ResidueIndex(const selection::SelectedModel& selected) {
        for (const auto& residue : selected.residues()) {
            residues_[key(residue.chain_name, residue.residue->seqid)].push_back(
                std::addressof(residue));
        }
    }

    [[nodiscard]] auto find(const std::string_view chain_name, const ::gemmi::SeqId& seqid) const
        -> std::span<const selection::SelectedResidue* const> {
        const auto found = residues_.find(key(chain_name, seqid));
        if (found == residues_.end()) {
            return {};
        }
        return found->second;
    }

  private:
    using Key = std::tuple<std::string_view, int, char>;

    // gemmi compares insertion codes ignoring case and treats a blank code as no code.
    [[nodiscard]] static auto key(const std::string_view chain_name, const ::gemmi::SeqId& seqid)
        -> Key {
        return {chain_name, seqid.num.value,
                static_cast<char>(static_cast<unsigned char>(seqid.icode) | 0x20U)};
    }

    std::map<Key, std::vector<const selection::SelectedResidue*>> residues_;
};

struct PartnerMatch {
    std::optional<std::size_t> atom;
    std::size_t residue_count = 0;
};

[[nodiscard]] auto match_partner(const ResidueIndex& residues,
                                 const selection::SelectedModel& selected,
                                 const ::gemmi::AtomAddress& address) -> PartnerMatch {
    PartnerMatch match;
    for (const auto* residue : residues.find(address.chain_name, address.res_id.seqid)) {
        if (!address.res_id.name.empty() && residue->component_name != address.res_id.name) {
            continue;
        }
        // The atom determines whether the partner is present, not which residue it belongs to.
        ++match.residue_count;
        const auto index = residue->component().find_atom(address.atom_name);
        if (index &&
            (address.altloc == '\0' || selected.atoms()[*index]->altloc_matches(address.altloc))) {
            match.atom = index;
        }
    }
    return match;
}

// Connection bond orders by _struct_conn.id; a connection without one is a single bond.
using ConnectionOrders = std::unordered_map<std::string, SourceOrder>;

// Adds covalent and disulfide connections between selected atoms and returns how many were skipped
// because their partners are in different symmetry images. The second partner of such a connection
// is an atom of a symmetry-related copy, which the molecule does not contain; its counterpart in
// the file is a different atom.
[[nodiscard]] auto add_structure_connections(BondAccumulator& bonds,
                                             const ::gemmi::Structure& structure,
                                             const selection::SelectedModel& selected,
                                             const ConnectionOrders& orders = {}) -> std::size_t {
    const auto is_bond = [](const ::gemmi::Connection& connection) {
        return connection.type == ::gemmi::Connection::Covale ||
               connection.type == ::gemmi::Connection::Disulf;
    };
    if (structure.models.empty() || !std::ranges::any_of(structure.connections, is_bond)) {
        return 0;
    }

    const auto residues = ResidueIndex{selected};
    auto skipped = std::size_t{0};
    for (const auto& connection : structure.connections) {
        if (!is_bond(connection)) {
            continue;
        }

        const auto first = match_partner(residues, selected, connection.partner1);
        const auto second = match_partner(residues, selected, connection.partner2);
        if (!first.atom || !second.atom) {
            continue;
        }
        if (connection.asu == ::gemmi::Asu::Different) {
            ++skipped;
            continue;
        }
        if (first.residue_count > 1 || second.residue_count > 1) {
            const auto& address =
                first.residue_count > 1 ? connection.partner1 : connection.partner2;
            const auto residue =
                detail::residue_description(address.res_id.name.empty() ? "?" : address.res_id.name,
                                            address.chain_name, address.res_id.seqid.str(), {});
            throw ParseError{ParseErrorKind::malformed,
                             "connection to atom " + address.atom_name +
                                 ": multiple selected residues match " + residue,
                             {}};
        }
        const auto order = orders.find(connection.name);
        bonds.add(*first.atom, *second.atom,
                  order == orders.end() ? SourceOrder{.order = core::BondOrder::SINGLE}
                                        : order->second);
    }
    return skipped;
}

// Whether two residues can be linked in a polymer of the given kind: both have templates of that
// kind, are in one chain, and have consecutive sequence numbers or insertion codes.
[[nodiscard]] auto sequentially_adjacent(const selection::SelectedResidue& previous,
                                         const selection::SelectedResidue& current,
                                         const component_templates::ComponentKind kind) -> bool {
    const auto* const previous_template = component_templates::find(previous.component().name);
    const auto* const current_template = component_templates::find(current.component().name);
    if (previous_template == nullptr || current_template == nullptr ||
        previous_template->kind != kind || current_template->kind != kind ||
        previous.chain_name != current.chain_name) {
        return false;
    }

    const auto& previous_id = previous.residue->seqid;
    const auto& current_id = current.residue->seqid;
    // Normal peptide/nucleotide progression increments the residue sequence number.
    if (previous_id.num.value + 1 == current_id.num.value) {
        return true;
    }
    // An insertion-code chain keeps the parent sequence number and progresses from the
    // unmarked base residue through consecutive codes (for example, 82, 82A, 82B).
    if (previous_id.num.value != current_id.num.value) {
        return false;
    }
    if (previous_id.icode == ' ') {
        return current_id.icode == 'A';
    }
    return current_id.icode == static_cast<char>(previous_id.icode + 1);
}

void add_sequential_bonds(BondAccumulator& bonds,
                          const std::span<const selection::SelectedResidue> residues,
                          const component_templates::ComponentKind kind,
                          const std::string_view previous_atom,
                          const std::string_view current_atom) {
    for (std::size_t index = 1; index < residues.size(); ++index) {
        const auto& previous = residues[index - 1];
        const auto& current = residues[index];
        if (!sequentially_adjacent(previous, current, kind)) {
            continue;
        }

        const auto first = previous.component().find_atom(previous_atom);
        const auto second = current.component().find_atom(current_atom);
        if (first.has_value() && second.has_value()) {
            bonds.add(*first, *second, core::BondOrder::SINGLE);
        }
    }
}

[[nodiscard]] auto bond_order(const std::string& value) -> std::optional<SourceOrder> {
    try {
        switch (::gemmi::bond_type_from_string(value)) {
        case ::gemmi::BondType::Single:
            return SourceOrder{.order = core::BondOrder::SINGLE};
        case ::gemmi::BondType::Double:
            return SourceOrder{.order = core::BondOrder::DOUBLE};
        case ::gemmi::BondType::Triple:
            return SourceOrder{.order = core::BondOrder::TRIPLE};
        case ::gemmi::BondType::Aromatic:
            return SourceOrder{.order = core::BondOrder::SINGLE, .aromatic = true};
        case ::gemmi::BondType::Deloc:
        case ::gemmi::BondType::Metal:
        case ::gemmi::BondType::Unspec:
            return std::nullopt;
        }
    } catch (const std::out_of_range&) {
        return std::nullopt;
    }

    return std::nullopt;
}

void add_component_bond(BondAccumulator& result, const detail::ComponentView component,
                        const auto& bond, const auto order) {
    const auto first = component.find_atom(bond.first);
    const auto second = component.find_atom(bond.second);
    if (first.has_value() && second.has_value()) {
        result.add(*first, *second, order);
    }
}

[[nodiscard]] auto
residue_protonation(const component_templates::ComponentTemplate& component_template,
                    const detail::ComponentView component, const selection::SelectedModel& model)
    -> std::span<const component_templates::BondTemplate> {
    return component_templates::protonation_bonds(
        component_template,
        [&](const std::string_view name) -> std::optional<component_templates::MotifAtom> {
            const auto index = component.find_atom(name);
            if (!index.has_value()) {
                return std::nullopt;
            }
            const auto& source = *model.atoms()[*index];
            return component_templates::MotifAtom{.atomic_number = source.element.atomic_number(),
                                                  .formal_charge = source.charge};
        });
}

using AtomPair = std::pair<std::size_t, std::size_t>;

[[nodiscard]] auto atom_pair(const core::Bond& bond) -> AtomPair {
    return {std::min(bond.first_atom_index(), bond.second_atom_index()),
            std::max(bond.first_atom_index(), bond.second_atom_index())};
}

// Whether a file may use a residue name for any compound. The wwPDB reserves these names for
// ligands that have no entry yet, so an entry of that name need not be the file's ligand.
[[nodiscard]] auto temporary_name(const std::string_view name) -> bool {
    const auto digit = [](const char character) { return character >= '0' && character <= '9'; };
    return name == "LIG" || name == "DRG" || name == "INH" ||
           (name.size() == 2 && digit(name[0]) && digit(name[1]));
}

// The dictionary entry of a residue, with the atoms of the entry that the residue has.
struct ResidueEntry {
    const DictionaryComponent* component;
    // For each atom of the entry, the selected atom of that name, if any.
    std::vector<std::optional<std::size_t>> atoms;
};

// The dictionary entries of residues. A structure repeats few component names in many residues, so
// each component is looked up once.
class DictionaryEntries {
  public:
    explicit DictionaryEntries(const ComponentDictionary& dictionary) : dictionary_{&dictionary} {}

    // The entry to use for a residue, if any. A residue with a temporary name is not looked up, and
    // one with an atom of another element than the entry's atom of that name does not use the
    // entry, which then describes other chemistry under the same names; both are recorded as
    // rejected. An element that the entry or the file does not give contradicts nothing, and atoms
    // of the entry may be absent from the residue.
    [[nodiscard]] auto find(const selection::SelectedResidue& residue,
                            const selection::SelectedModel& model) -> std::optional<ResidueEntry> {
        const auto component = residue.component();
        if (temporary_name(component.name)) {
            reject(residue, "a temporary ligand name");
            return std::nullopt;
        }
        auto found = entries_.find(component.name);
        if (found == entries_.end()) {
            found = entries_.emplace(component.name, dictionary_->find(component.name)).first;
        }
        if (!found->second.has_value()) {
            return std::nullopt;
        }
        auto entry = ResidueEntry{.component = std::addressof(*found->second), .atoms = {}};
        for (const auto& atom : entry.component->atoms) {
            const auto index = entry.atoms.emplace_back(component.find_atom(atom.name));
            if (index.has_value()) {
                const auto atomic_number = model.atoms()[*index]->element.atomic_number();
                if (atom.atomic_number != 0 && atomic_number != 0 &&
                    atom.atomic_number != atomic_number) {
                    reject(residue, "another element for atom " + atom.name + " in its entry");
                    return std::nullopt;
                }
            }
        }
        return entry;
    }

    // Records that a residue takes nothing from the dictionary, although it has no other source.
    void reject(const selection::SelectedResidue& residue, const std::string_view reason) {
        constexpr std::size_t example_limit = 3;
        if (++rejected_count_ <= example_limit) {
            rejected_examples_ += rejected_count_ == 1 ? "" : ", ";
            rejected_examples_ +=
                detail::residue_description(residue.component_name, residue.chain_name,
                                            residue.residue->seqid.str(), {}) +
                " (" + std::string{reason} + ")";
        }
    }

    // The record warning that names the rejected residues, if any.
    [[nodiscard]] auto issue() const -> std::optional<ImportIssue> {
        if (rejected_count_ == 0) {
            return std::nullopt;
        }
        return ImportIssue{.kind = ImportIssueKind::component_entries_rejected,
                           .message = "the component dictionary " + dictionary_->version() +
                                      " was not used for " + std::to_string(rejected_count_) +
                                      (rejected_count_ == 1 ? " residue" : " residues") +
                                      ", for example " + rejected_examples_,
                           .line = std::nullopt};
    }

  private:
    std::size_t rejected_count_ = 0;
    std::string rejected_examples_;
    const ComponentDictionary* dictionary_;
    std::unordered_map<std::string_view, std::optional<DictionaryComponent>> entries_;
};

// The bonds of an entry between the atoms present in its residue.
[[nodiscard]] auto entry_bonds(const ResidueEntry& entry) -> std::vector<core::Bond> {
    auto result = std::vector<core::Bond>{};
    for (const auto& bond : entry.component->bonds) {
        const auto first = entry.atoms[bond.first];
        const auto second = entry.atoms[bond.second];
        if (first.has_value() && second.has_value() && *first != *second) {
            result.emplace_back(*first, *second, bond.order);
        }
    }
    return result;
}

// How an entry differs from the bonds that the file states within a residue, or nothing when it
// describes them, so that the file's bonds can take the entry's orders. Equal names alone do not
// show that: an entry renamed since the file was written, or the entry of another compound with the
// same residue name, shares some. The entry must name every non-hydrogen atom of the residue and
// bond these atoms in the same pairs as the file, with no pair more or less. A hydrogen it names
// must have the partner the file gives it; one it does not name, such as one a preparation tool
// named, is left out. A match does not prove that the residue is the entry's compound, or in the
// entry's protonation state.
[[nodiscard]] auto entry_mismatch(const selection::SelectedResidue& residue,
                                  const selection::SelectedModel& model, const ResidueEntry& entry,
                                  const std::map<AtomPair, core::BondOrder>& entry_orders,
                                  const std::span<const AtomPair> file_pairs)
    -> std::optional<std::string> {
    auto named = std::unordered_set<std::size_t>{};
    for (const auto& atom : entry.atoms) {
        if (atom.has_value()) {
            named.insert(*atom);
        }
    }
    const auto is_hydrogen = [&](const std::size_t atom) {
        return model.atoms()[atom]->element.is_hydrogen();
    };
    const auto heavy = [&](const AtomPair& pair) {
        return !is_hydrogen(pair.first) && !is_hydrogen(pair.second);
    };
    const auto unnamed = std::ranges::find_if(residue.atoms, [&](const auto& atom) {
        return !is_hydrogen(atom.index) && !named.contains(atom.index);
    });
    if (unnamed != residue.atoms.end()) {
        return "no atom " + unnamed->name + " in its entry";
    }
    // Every pair of the file is one of the entry's, or bonds a hydrogen that the entry does not
    // name. Such a hydrogen is not a partner for a hydrogen the entry names, which has its own.
    const auto named_hydrogen = [&](const std::size_t atom) {
        return is_hydrogen(atom) && named.contains(atom);
    };
    const auto file_pairs_fit = std::ranges::all_of(file_pairs, [&](const AtomPair& pair) {
        return entry_orders.contains(pair) ||
               ((!named.contains(pair.first) || !named.contains(pair.second)) &&
                !named_hydrogen(pair.first) && !named_hydrogen(pair.second));
    });
    // The file's pairs are distinct, so equal counts leave no pair of the entry without one.
    if (!file_pairs_fit || std::ranges::count_if(file_pairs, heavy) !=
                               std::ranges::count_if(entry_orders, heavy,
                                                     [](const auto& bond) { return bond.first; })) {
        return "other bonds in its entry than in the file";
    }
    return std::nullopt;
}

struct DictionaryBonds {
    std::vector<core::Bond> added;
    std::optional<ImportIssue> issue;
};

// What the residues without a built-in template take from the dictionary. Its atom names need not
// mean the atoms of a residue that the file bonds itself, so there are three cases:
// - a residue that the file states no bond within takes the bonds of its entry, which are returned;
// - a residue that the file bonds takes orders alone, for its bonds marked in order_unknown, and
//   only when its entry matches it: these bonds of `file` get the orders and lose their marks;
// - any other residue takes nothing.
// A _chem_comp_bond row of the residue's component counts as a bond of the file, also when no bond
// was read from it; a connection to another residue does not.
[[nodiscard]] auto apply_dictionary(const selection::SelectedModel& model,
                                    const ComponentDictionary& dictionary, ExplicitBonds& file,
                                    std::vector<bool>& order_unknown) -> DictionaryBonds {
    // The file's bonds within each residue, by their indices in file.bonds.
    const auto& sites = model.sites();
    auto file_bonds = std::unordered_map<const ::gemmi::Residue*, std::vector<std::size_t>>{};
    for (std::size_t index = 0; index < file.bonds.size(); ++index) {
        const auto* const residue = sites[file.bonds[index].first_atom_index()].residue;
        if (residue == sites[file.bonds[index].second_atom_index()].residue) {
            file_bonds[residue].push_back(index);
        }
    }

    auto entries = DictionaryEntries{dictionary};
    BondAccumulator added;
    for (const auto& residue : model.residues()) {
        if (component_templates::find(residue.component_name) != nullptr) {
            continue;
        }
        const auto bonded = file_bonds.find(residue.residue);
        if (bonded == file_bonds.end()) {
            if (!file.components_with_bond_rows.contains(residue.component_name)) {
                if (const auto entry = entries.find(residue, model)) {
                    for (const auto& bond : entry_bonds(*entry)) {
                        added.add(bond.first_atom_index(), bond.second_atom_index(), bond.order());
                    }
                }
            }
            continue;
        }
        const auto unknown = [&](const std::size_t index) { return order_unknown[index]; };
        if (std::ranges::none_of(bonded->second, unknown)) {
            continue;
        }
        const auto entry = entries.find(residue, model);
        if (!entry.has_value()) {
            continue;
        }
        auto entry_orders = std::map<AtomPair, core::BondOrder>{};
        for (const auto& bond : entry_bonds(*entry)) {
            entry_orders.emplace(atom_pair(bond), bond.order());
        }
        auto file_pairs = std::vector<AtomPair>{};
        for (const auto index : bonded->second) {
            file_pairs.push_back(atom_pair(file.bonds[index]));
        }
        if (const auto mismatch =
                entry_mismatch(residue, model, *entry, entry_orders, file_pairs)) {
            entries.reject(residue, *mismatch);
            continue;
        }
        for (const auto index : bonded->second | std::views::filter(unknown)) {
            auto& bond = file.bonds[index];
            bond = core::Bond{bond.first_atom_index(), bond.second_atom_index(),
                              entry_orders.at(atom_pair(bond))};
            order_unknown[index] = false;
        }
    }
    return {.added = added.take(), .issue = entries.issue()};
}

[[nodiscard]] auto assign_template_bonds(const selection::SelectedModel& model)
    -> std::vector<core::Bond> {
    BondAccumulator result;
    const auto& residues = model.residues();

    for (std::size_t residue_index = 0; residue_index < residues.size(); ++residue_index) {
        const auto& residue = residues[residue_index];
        const auto component = residue.component();
        const auto* const component_template = component_templates::find(component.name);
        if (component_template == nullptr) {
            continue;
        }

        // These orders go in first, so the accumulator ignores the template's for the same pairs.
        for (const auto& bond : residue_protonation(*component_template, component, model)) {
            add_component_bond(result, component, bond, bond.order);
        }
        for (const auto& bond : component_template->bonds) {
            add_component_bond(result, component, bond, bond.order);
        }
        // Only an amino acid that the templates do not link to a previous residue is an N-terminus,
        // as add_sequential_bonds() decides below. Hydride also names H1 a hydrogen on the carboxyl
        // carbon of a C-terminal residue without OXT.
        const auto peptide_linked =
            residue_index > 0 &&
            sequentially_adjacent(residues[residue_index - 1], residue,
                                  component_templates::ComponentKind::amino_acid) &&
            residues[residue_index - 1].component().find_atom("C").has_value();
        const auto n_terminal =
            component_template->kind == component_templates::ComponentKind::amino_acid &&
            !peptide_linked;
        if (n_terminal) {
            if (const auto nitrogen = component.find_atom("N")) {
                // Supplement the base template with the remaining terminal hydrogen names. Hydride
                // names the second hydrogen of an N-terminal proline H1 rather than H2 or H3.
                static constexpr auto proline_names =
                    std::array<std::string_view, 3>{"H1", "H2", "H3"};
                static constexpr auto other_names = std::array<std::string_view, 2>{"H1", "H3"};
                const auto names = component.name == "PRO"
                                       ? std::span<const std::string_view>{proline_names}
                                       : std::span<const std::string_view>{other_names};
                for (const auto name : names) {
                    if (const auto hydrogen = component.find_atom(name)) {
                        result.add(*nitrogen, *hydrogen, core::BondOrder::SINGLE);
                    }
                }
            }
        }
    }

    add_sequential_bonds(result, residues, component_templates::ComponentKind::amino_acid, "C",
                         "N");
    add_sequential_bonds(result, residues, component_templates::ComponentKind::nucleotide, "O3'",
                         "P");

    return result.take();
}

} // namespace

auto explicit_pdb(const ::gemmi::Structure& structure, const selection::SelectedModel& model)
    -> ExplicitBonds {
    if (structure.models.empty()) {
        return {};
    }

    BondAccumulator result;
    const auto skipped = add_structure_connections(result, structure, model);
    auto unknown = std::vector<std::size_t>{};

    for (const auto& [serial, partners] : structure.conect_map) {
        const auto first = model.atom_index_by_serial(serial);
        if (!first.has_value()) {
            continue;
        }

        for (const auto partner : partners) {
            const auto second = model.atom_index_by_serial(partner);
            if (second.has_value()) {
                const auto index = result.add(*first, *second, core::BondOrder::SINGLE);
                const auto& first_site = model.sites()[*first];
                const auto& second_site = model.sites()[*second];
                // A bond to hydrogen can only be single, so its order is known.
                if (index && first_site.residue == second_site.residue &&
                    !first_site.atom->element.is_hydrogen() &&
                    !second_site.atom->element.is_hydrogen()) {
                    unknown.push_back(*index);
                }
            }
        }
    }

    return {.bonds = result.take(),
            .orders = ExplicitBondOrders::connectivity_only,
            .skipped_symmetry_connections = skipped,
            .unknown_order_indices = std::move(unknown)};
}

auto explicit_mmcif(const ::gemmi::Structure& structure, ::gemmi::cif::Block& block,
                    const selection::SelectedModel& model, const BondStrategy strategy)
    -> ExplicitBonds {
    if (structure.models.empty()) {
        return {};
    }

    const auto& residues = model.residues();
    BondAccumulator result;

    struct ParsedComponentBond {
        std::string first;
        std::string second;
        SourceOrder order;
    };

    std::unordered_map<std::string, std::vector<ParsedComponentBond>> component_bonds;
    // A component is listed here once the file has a row for it, also when no row gives a bond.
    auto rows =
        block.find("_chem_comp_bond.", {"comp_id", "atom_id_1", "atom_id_2", "?value_order"});
    for (const auto row : rows) {
        auto& bonds = component_bonds[row.str(0)];
        if (!row.has(3)) {
            continue;
        }
        if (const auto order = bond_order(row.str(3))) {
            bonds.emplace_back(
                ParsedComponentBond{.first = row.str(1), .second = row.str(2), .order = *order});
        }
    }

    for (const auto& residue : residues) {
        const auto component = residue.component();
        const auto found = component_bonds.find(std::string{component.name});
        if (found == component_bonds.end()) {
            continue;
        }
        // The rows describe a component, so they cannot follow the protonation of one residue.
        // Under hybrid, a bond that the residue's template has takes the template's order.
        const auto* const component_template =
            strategy == BondStrategy::hybrid ? component_templates::find(component.name) : nullptr;
        const auto protonation = component_template != nullptr
                                     ? residue_protonation(*component_template, component, model)
                                     : std::span<const component_templates::BondTemplate>{};
        for (const auto& bond : found->second) {
            auto order = bond.order;
            if (component_template != nullptr) {
                if (const auto known = component_templates::bond_order(
                        *component_template, protonation, bond.first, bond.second)) {
                    order = SourceOrder{.order = *known};
                }
            }
            add_component_bond(result, component, bond, order);
        }
    }

    // Gemmi keeps each connection's _struct_conn.id but not its bond order.
    auto connection_orders = ConnectionOrders{};
    for (const auto row : block.find("_struct_conn.", {"id", "?pdbx_value_order"})) {
        if (!row.has2(1)) {
            continue;
        }
        if (const auto order = bond_order(row.str(1))) {
            connection_orders.emplace(row.str(0), *order);
        }
    }
    const auto skipped = add_structure_connections(result, structure, model, connection_orders);

    auto components_with_bond_rows = std::unordered_set<std::string>{};
    for (const auto& component : component_bonds | std::views::keys) {
        components_with_bond_rows.insert(component);
    }
    return {.bonds = result.take(),
            .skipped_symmetry_connections = skipped,
            .aromatic_bonds = result.aromatic_bonds(),
            .components_with_bond_rows = std::move(components_with_bond_rows)};
}

auto symmetry_connections_issue(const std::size_t count) -> ImportIssue {
    return ImportIssue{.kind = ImportIssueKind::symmetry_connections_ignored,
                       .message = std::to_string(count) +
                                  (count == 1 ? " bond to an atom of a symmetry-related copy was"
                                              : " bonds to atoms of symmetry-related copies were") +
                                  " ignored; only the atoms in the file are calculated",
                       .line = std::nullopt};
}

auto assign(const selection::SelectedModel& model, const BondStrategy strategy,
            const ComponentDictionary* const dictionary, ExplicitBonds explicit_bonds)
    -> Assignment {
    if (explicit_bonds.orders == ExplicitBondOrders::from_source &&
        !explicit_bonds.unknown_order_indices.empty()) {
        throw std::logic_error{"unknown-order indices require connectivity-only input"};
    }
    switch (strategy) {
    case BondStrategy::none:
        return {};
    case BondStrategy::templates: {
        auto result = Assignment{.bonds = assign_template_bonds(model)};
        if (dictionary != nullptr) {
            // The file's bonds are not read, so every residue without a built-in template takes
            // the bonds of its entry. The residues are other ones, so no pair is given twice.
            auto no_bonds = ExplicitBonds{};
            auto no_orders = std::vector<bool>{};
            auto from_dictionary = apply_dictionary(model, *dictionary, no_bonds, no_orders);
            std::ranges::move(from_dictionary.added, std::back_inserter(result.bonds));
            result.dictionary_issue = std::move(from_dictionary.issue);
        }
        return result;
    }
    case BondStrategy::explicit_bonds:
        return {.bonds = std::move(explicit_bonds.bonds),
                .unknown_order_indices = std::move(explicit_bonds.unknown_order_indices)};
    case BondStrategy::hybrid: {
        // The first source to define a pair sets its order. Bonds with source orders (mmCIF) come
        // first and keep their place in the record; explicit_mmcif() has given the component bonds
        // that a template has the template's orders. Connectivity-only records (PDB CONECT) only
        // add pairs templates do not cover. The dictionary's bonds come last, and the orders it
        // gives are put on the file's bonds before the sources are merged.
        auto result = Assignment{};
        auto order_unknown = std::vector<bool>(explicit_bonds.bonds.size());
        for (const auto index : explicit_bonds.unknown_order_indices) {
            order_unknown[index] = true;
        }
        auto from_dictionary = DictionaryBonds{};
        if (dictionary != nullptr) {
            from_dictionary = apply_dictionary(model, *dictionary, explicit_bonds, order_unknown);
            result.dictionary_issue = std::move(from_dictionary.issue);
        }
        const auto file_first = explicit_bonds.orders == ExplicitBondOrders::from_source;
        auto first = assign_template_bonds(model);
        auto second = std::move(explicit_bonds.bonds);
        if (file_first) {
            std::swap(first, second);
        }

        BondAccumulator bonds;
        for (const auto& bond : first) {
            bonds.add(bond.first_atom_index(), bond.second_atom_index(), bond.order());
        }
        for (std::size_t i = 0; i < second.size(); ++i) {
            const auto& bond = second[i];
            const auto index =
                bonds.add(bond.first_atom_index(), bond.second_atom_index(), bond.order());
            // A bond that a template covers is not added again, and its order is the template's.
            if (!file_first && order_unknown[i] && index.has_value()) {
                result.unknown_order_indices.push_back(*index);
            }
        }
        for (const auto& bond : from_dictionary.added) {
            bonds.add(bond.first_atom_index(), bond.second_atom_index(), bond.order());
        }
        result.bonds = bonds.take();
        return result;
    }
    }

    std::abort();
}

} // namespace chargefw::io::gemmi::bonds
