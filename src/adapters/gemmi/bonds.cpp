#include "bonds.h"

#include "component_templates.h"
#include "selection.h"

#include <gemmi/chemcomp.hpp>
#include <gemmi/cif.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdlib>
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

namespace chargefw::adapters::gemmi::bonds {
namespace {

// A bond order read from a file. Aromatic bonds are read as single bonds.
struct SourceOrder {
    core::BondOrder order;
    bool aromatic = false;
};

class BondAccumulator {
  public:
    void add(const std::size_t first, const std::size_t second, const core::BondOrder order) {
        add(first, second, SourceOrder{.order = order});
    }

    // The first order given for a pair is kept.
    void add(const std::size_t first, const std::size_t second, const SourceOrder order) {
        if (first == second) {
            return;
        }

        const auto key =
            BondKey{.first = std::min(first, second), .second = std::max(first, second)};
        if (seen_.insert(key).second) {
            bonds_.emplace_back(first, second, order.order);
            if (order.aromatic) {
                aromatic_bonds_.push_back(bonds_.size() - 1);
            }
        }
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

// Residues of a model by chain name and sequence ID, in model order, so that each connection
// partner is found without scanning every chain and residue.
class ResidueIndex {
  public:
    explicit ResidueIndex(const ::gemmi::Model& model) {
        for (const auto& chain : model.chains) {
            for (const auto& residue : chain.residues) {
                residues_[key(chain.name, residue.seqid)].push_back(std::addressof(residue));
            }
        }
    }

    [[nodiscard]] auto find(const std::string_view chain_name, const ::gemmi::SeqId& seqid) const
        -> std::span<const ::gemmi::Residue* const> {
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

    std::map<Key, std::vector<const ::gemmi::Residue*>> residues_;
};

[[nodiscard]] auto selected_atom_index(const ResidueIndex& residues,
                                       const selection::SelectedModel& selected,
                                       const ::gemmi::AtomAddress& address)
    -> std::optional<std::size_t> {
    const auto candidates = residues.find(address.chain_name, address.res_id.seqid);
    // As Model::find_cra(address, true): the first residue with the given name decides.
    for (const auto* residue : candidates) {
        if (residue->name != address.res_id.name) {
            continue;
        }
        if (!address.atom_name.empty()) {
            if (const auto* atom = residue->find_atom(address.atom_name, address.altloc)) {
                return selected.atom_index(atom);
            }
        }
        break;
    }
    if (!address.res_id.name.empty() && address.altloc != '\0') {
        return std::nullopt;
    }

    for (const auto* residue : candidates) {
        if (!address.res_id.name.empty() && residue->name != address.res_id.name) {
            continue;
        }
        for (const auto& atom : residue->atoms) {
            if (atom.name == address.atom_name &&
                (address.altloc == '\0' || atom.altloc == address.altloc)) {
                return selected.atom_index(std::addressof(atom));
            }
        }
    }
    return std::nullopt;
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

    const auto residues = ResidueIndex{structure.models.front()};
    auto skipped = std::size_t{0};
    for (const auto& connection : structure.connections) {
        if (!is_bond(connection)) {
            continue;
        }

        const auto first = selected_atom_index(residues, selected, connection.partner1);
        const auto second = selected_atom_index(residues, selected, connection.partner2);
        if (!first.has_value() || !second.has_value()) {
            continue;
        }
        if (connection.asu == ::gemmi::Asu::Different) {
            ++skipped;
            continue;
        }
        const auto order = orders.find(connection.name);
        bonds.add(*first, *second,
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

void add_component_bonds(BondAccumulator& result, const detail::ComponentView component,
                         const auto& bonds) {
    for (const auto& bond : bonds) {
        const auto first = component.find_atom(bond.first);
        const auto second = component.find_atom(bond.second);
        if (first.has_value() && second.has_value()) {
            result.add(*first, *second, bond.order);
        }
    }
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
        const auto atom =
            [&](const std::string_view name) -> std::optional<component_templates::MotifAtom> {
            const auto index = component.find_atom(name);
            if (!index.has_value()) {
                return std::nullopt;
            }
            const auto& source = *model.atoms()[*index];
            return component_templates::MotifAtom{.atomic_number = source.element.atomic_number(),
                                                  .formal_charge = source.charge};
        };
        add_component_bonds(result, component,
                            component_templates::protonation_bonds(*component_template, atom));
        add_component_bonds(result, component, component_template->bonds);
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

    for (const auto& [serial, partners] : structure.conect_map) {
        const auto first = model.atom_index_by_serial(serial);
        if (!first.has_value()) {
            continue;
        }

        for (const auto partner : partners) {
            const auto second = model.atom_index_by_serial(partner);
            if (second.has_value()) {
                result.add(*first, *second, core::BondOrder::SINGLE);
            }
        }
    }

    return {.bonds = result.take(), .skipped_symmetry_connections = skipped};
}

auto explicit_mmcif(const ::gemmi::Structure& structure, ::gemmi::cif::Block& block,
                    const selection::SelectedModel& model) -> ExplicitBonds {
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
    auto rows =
        block.find("_chem_comp_bond.", {"comp_id", "atom_id_1", "atom_id_2", "value_order"});
    for (const auto row : rows) {
        auto component = row.str(0);
        auto first = row.str(1);
        auto second = row.str(2);
        const auto order_value = row.str(3);
        if (const auto order = bond_order(order_value)) {
            component_bonds[std::move(component)].emplace_back(ParsedComponentBond{
                .first = std::move(first), .second = std::move(second), .order = *order});
        }
    }

    for (const auto& residue : residues) {
        const auto component = residue.component();
        const auto found = component_bonds.find(std::string{component.name});
        if (found != component_bonds.end()) {
            add_component_bonds(result, component, found->second);
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

    return {.bonds = result.take(),
            .skipped_symmetry_connections = skipped,
            .aromatic_bonds = result.aromatic_bonds()};
}

auto symmetry_connections_issue(const std::size_t count) -> ImportIssue {
    return ImportIssue{.kind = ImportIssueKind::symmetry_connections_ignored,
                       .message = std::to_string(count) +
                                  (count == 1 ? " bond to an atom of a symmetry-related copy was"
                                              : " bonds to atoms of symmetry-related copies were") +
                                  " ignored; only the atoms in the file are calculated.",
                       .line = std::nullopt};
}

auto assign(const selection::SelectedModel& model, const BondStrategy strategy,
            std::vector<core::Bond> explicit_bonds, const ExplicitBondOrders explicit_orders)
    -> std::vector<core::Bond> {
    switch (strategy) {
    case BondStrategy::none:
        return {};
    case BondStrategy::templates:
        return assign_template_bonds(model);
    case BondStrategy::explicit_bonds:
        return explicit_bonds;
    case BondStrategy::hybrid: {
        // The first source to define a pair sets its order. Source orders take precedence over
        // templates; connectivity-only records (PDB CONECT) only add pairs templates do not cover.
        auto first = assign_template_bonds(model);
        auto second = std::move(explicit_bonds);
        if (explicit_orders == ExplicitBondOrders::from_source) {
            std::swap(first, second);
        }

        BondAccumulator bonds;
        for (const auto& bond : first) {
            bonds.add(bond.first_atom_index(), bond.second_atom_index(), bond.order());
        }
        for (const auto& bond : second) {
            bonds.add(bond.first_atom_index(), bond.second_atom_index(), bond.order());
        }

        return bonds.take();
    }
    }

    std::abort();
}

} // namespace chargefw::adapters::gemmi::bonds
