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
#include <functional>
#include <iterator>
#include <map>
#include <memory>
#include <optional>
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

// Which selected residues the file states a bond within.
class InternalBonds {
  public:
    explicit InternalBonds(const selection::SelectedModel& model)
        : model_{&model}, stated_(model.residues().size()) {}

    // Notes a pair of atoms the file bonds, whether or not the bond is new.
    void note(const std::size_t first, const std::size_t second) {
        const auto* const residue = model_->sites()[first].residue;
        if (first != second && residue == model_->sites()[second].residue) {
            noted_.insert(residue);
        }
    }

    // Notes a residue by its index in the model.
    void note_residue(const std::size_t residue) {
        stated_[residue] = true;
    }

    [[nodiscard]] auto take() -> std::vector<bool> {
        const auto& residues = model_->residues();
        for (std::size_t index = 0; index < residues.size() && !noted_.empty(); ++index) {
            if (noted_.contains(residues[index].residue)) {
                stated_[index] = true;
            }
        }
        return std::move(stated_);
    }

  private:
    const selection::SelectedModel* model_;
    std::vector<bool> stated_;
    // The residues of the noted pairs. Files state few bonds this way, so the residues are looked
    // up once at the end, not indexed by atom beforehand.
    std::unordered_set<const ::gemmi::Residue*> noted_;
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
[[nodiscard]] auto add_structure_connections(BondAccumulator& bonds, InternalBonds& internal,
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
        internal.note(*first.atom, *second.atom);
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

// The dictionary's bonds between the atoms present in the residues that have no built-in template
// and that use_residue accepts by their index in the model.
[[nodiscard]] auto dictionary_bonds(const selection::SelectedModel& model,
                                    const ComponentDictionary& dictionary,
                                    const std::function<bool(std::size_t)>& use_residue)
    -> std::vector<core::Bond> {
    BondAccumulator result;
    // A structure repeats few component names in many residues.
    auto entries = std::unordered_map<std::string_view, std::optional<DictionaryComponent>>{};
    const auto& residues = model.residues();
    for (std::size_t index = 0; index < residues.size(); ++index) {
        const auto component = residues[index].component();
        if (component_templates::find(component.name) != nullptr || !use_residue(index)) {
            continue;
        }
        auto entry = entries.find(component.name);
        if (entry == entries.end()) {
            entry = entries.emplace(component.name, dictionary.find(component.name)).first;
        }
        if (!entry->second.has_value()) {
            continue;
        }
        const auto& atoms = entry->second->atoms;
        for (const auto& bond : entry->second->bonds) {
            const auto first = component.find_atom(atoms[bond.first].name);
            const auto second = component.find_atom(atoms[bond.second].name);
            if (first.has_value() && second.has_value()) {
                result.add(*first, *second, bond.order);
            }
        }
    }
    return result.take();
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
    auto internal = InternalBonds{model};
    const auto skipped = add_structure_connections(result, internal, structure, model);
    auto unknown = std::vector<std::size_t>{};

    for (const auto& [serial, partners] : structure.conect_map) {
        const auto first = model.atom_index_by_serial(serial);
        if (!first.has_value()) {
            continue;
        }

        for (const auto partner : partners) {
            const auto second = model.atom_index_by_serial(partner);
            if (second.has_value()) {
                internal.note(*first, *second);
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
            .skipped_symmetry_connections = skipped,
            .unknown_order_indices = std::move(unknown),
            .residues_with_file_bonds = internal.take()};
}

auto explicit_mmcif(const ::gemmi::Structure& structure, ::gemmi::cif::Block& block,
                    const selection::SelectedModel& model, const BondStrategy strategy)
    -> ExplicitBonds {
    if (structure.models.empty()) {
        return {};
    }

    const auto& residues = model.residues();
    BondAccumulator result;
    auto internal = InternalBonds{model};

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

    for (std::size_t residue_index = 0; residue_index < residues.size(); ++residue_index) {
        const auto& residue = residues[residue_index];
        const auto component = residue.component();
        const auto found = component_bonds.find(std::string{component.name});
        if (found == component_bonds.end()) {
            continue;
        }
        internal.note_residue(residue_index);
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
    const auto skipped =
        add_structure_connections(result, internal, structure, model, connection_orders);

    return {.bonds = result.take(),
            .skipped_symmetry_connections = skipped,
            .aromatic_bonds = result.aromatic_bonds(),
            .residues_with_file_bonds = internal.take()};
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
            const ComponentDictionary* const dictionary, std::vector<core::Bond> explicit_bonds,
            const std::vector<bool>& residues_with_file_bonds,
            const ExplicitBondOrders explicit_orders,
            std::vector<std::size_t>& unknown_order_indices) -> std::vector<core::Bond> {
    if (explicit_orders == ExplicitBondOrders::from_source && !unknown_order_indices.empty()) {
        throw std::logic_error{"unknown-order indices require connectivity-only input"};
    }
    switch (strategy) {
    case BondStrategy::none:
        unknown_order_indices.clear();
        return {};
    case BondStrategy::templates: {
        unknown_order_indices.clear();
        auto bonds = assign_template_bonds(model);
        if (dictionary != nullptr) {
            // The residues are other ones, so no pair is given twice.
            std::ranges::move(
                dictionary_bonds(model, *dictionary, [](std::size_t) { return true; }),
                std::back_inserter(bonds));
        }
        return bonds;
    }
    case BondStrategy::explicit_bonds:
        return explicit_bonds;
    case BondStrategy::hybrid: {
        // The first source to define a pair sets its order. Bonds with source orders (mmCIF) come
        // first and keep their place in the record; explicit_mmcif() has given the component bonds
        // that a template has the template's orders. Connectivity-only records (PDB CONECT) only
        // add pairs templates do not cover. The dictionary comes last and bonds only the residues
        // that have no built-in template and no bond within them in the file: its atom names need
        // not mean the atoms of a residue the file bonds itself. In such a residue it gives only
        // the order of a pair whose order the file does not state.
        const auto file_bonded = [&](const std::size_t residue) {
            return residue < residues_with_file_bonds.size() && residues_with_file_bonds[residue];
        };
        auto from_dictionary = std::vector<core::Bond>{};
        auto dictionary_orders = std::map<std::pair<std::size_t, std::size_t>, core::BondOrder>{};
        if (dictionary != nullptr) {
            from_dictionary = dictionary_bonds(model, *dictionary, [&](const std::size_t residue) {
                return !file_bonded(residue);
            });
            if (!unknown_order_indices.empty()) {
                for (const auto& bond : dictionary_bonds(model, *dictionary, file_bonded)) {
                    dictionary_orders.emplace(
                        std::minmax(bond.first_atom_index(), bond.second_atom_index()),
                        bond.order());
                }
            }
        }
        auto first = assign_template_bonds(model);
        auto second = std::move(explicit_bonds);
        if (explicit_orders == ExplicitBondOrders::from_source) {
            std::swap(first, second);
        }

        BondAccumulator bonds;
        for (const auto& bond : first) {
            bonds.add(bond.first_atom_index(), bond.second_atom_index(), bond.order());
        }
        auto retained_unknown = std::vector<std::size_t>{};
        auto unknown = unknown_order_indices.begin();
        for (std::size_t i = 0; i < second.size(); ++i) {
            const auto& bond = second[i];
            const auto order_unknown = unknown != unknown_order_indices.end() && *unknown == i;
            const auto known = order_unknown
                                   ? dictionary_orders.find(std::minmax(bond.first_atom_index(),
                                                                        bond.second_atom_index()))
                                   : dictionary_orders.end();
            const auto index =
                bonds.add(bond.first_atom_index(), bond.second_atom_index(),
                          known != dictionary_orders.end() ? known->second : bond.order());
            if (order_unknown) {
                if (index && known == dictionary_orders.end()) {
                    retained_unknown.push_back(*index);
                }
                ++unknown;
            }
        }
        unknown_order_indices = std::move(retained_unknown);
        for (const auto& bond : from_dictionary) {
            bonds.add(bond.first_atom_index(), bond.second_atom_index(), bond.order());
        }

        return bonds.take();
    }
    }

    std::abort();
}

} // namespace chargefw::io::gemmi::bonds
