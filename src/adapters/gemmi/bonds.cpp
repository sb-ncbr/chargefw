#include "bonds.h"

#include "component_templates.h"
#include "selection.h"

#include <gemmi/chemcomp.hpp>
#include <gemmi/cif.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace chargefw::adapters::gemmi::bonds {
namespace {

class BondAccumulator {
  public:
    void add(const std::size_t first, const std::size_t second, const core::BondOrder order) {
        if (first == second) {
            return;
        }

        const auto key =
            BondKey{.first = std::min(first, second), .second = std::max(first, second)};
        if (seen_.insert(key).second) {
            bonds_.emplace_back(first, second, order);
        }
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
};

[[nodiscard]] auto selected_atom_index(const ::gemmi::Model& model,
                                       const selection::SelectedModel& selected,
                                       const ::gemmi::AtomAddress& address)
    -> std::optional<std::size_t> {
    const auto cra = model.find_cra(address, true);
    if (cra.atom != nullptr) {
        return selected.atom_index(cra.atom);
    }
    if (!address.res_id.name.empty() && address.altloc != '\0') {
        return std::nullopt;
    }

    for (const auto& chain : model.chains) {
        if (chain.name != address.chain_name) {
            continue;
        }
        for (const auto& residue : chain.residues) {
            if (residue.seqid != address.res_id.seqid ||
                (!address.res_id.name.empty() && residue.name != address.res_id.name)) {
                continue;
            }
            for (const auto& atom : residue.atoms) {
                if (atom.name == address.atom_name &&
                    (address.altloc == '\0' || atom.altloc == address.altloc)) {
                    return selected.atom_index(std::addressof(atom));
                }
            }
        }
    }
    return std::nullopt;
}

// Connection bond orders by _struct_conn.id; a connection without one is a single bond.
using ConnectionOrders = std::unordered_map<std::string, core::BondOrder>;

void add_structure_connections(BondAccumulator& bonds, const ::gemmi::Structure& structure,
                               const selection::SelectedModel& selected,
                               const ConnectionOrders& orders = {}) {
    if (structure.models.empty()) {
        return;
    }

    const auto& model = structure.models.front();
    for (const auto& connection : structure.connections) {
        if (connection.type != ::gemmi::Connection::Covale &&
            connection.type != ::gemmi::Connection::Disulf) {
            continue;
        }

        const auto first = selected_atom_index(model, selected, connection.partner1);
        const auto second = selected_atom_index(model, selected, connection.partner2);
        if (first.has_value() && second.has_value()) {
            const auto order = orders.find(connection.name);
            bonds.add(*first, *second,
                      order == orders.end() ? core::BondOrder::SINGLE : order->second);
        }
    }
}

void add_sequential_bonds(BondAccumulator& bonds,
                          const std::span<const selection::SelectedResidue> residues,
                          const component_templates::ComponentKind kind,
                          const std::string_view previous_atom,
                          const std::string_view current_atom) {
    const auto sequentially_adjacent = [](const ::gemmi::Residue& previous,
                                          const ::gemmi::Residue& current) -> bool {
        // Normal peptide/nucleotide progression increments the residue sequence number.
        if (previous.seqid.num.value + 1 == current.seqid.num.value) {
            return true;
        }
        // An insertion-code chain keeps the parent sequence number and progresses from the
        // unmarked base residue through consecutive codes (for example, 82, 82A, 82B).
        if (previous.seqid.num.value != current.seqid.num.value) {
            return false;
        }
        if (previous.seqid.icode == ' ') {
            return current.seqid.icode == 'A';
        }
        return current.seqid.icode == static_cast<char>(previous.seqid.icode + 1);
    };

    for (std::size_t index = 1; index < residues.size(); ++index) {
        const auto& previous = residues[index - 1];
        const auto& current = residues[index];
        const auto previous_component = previous.component();
        const auto current_component = current.component();
        const auto* const previous_template = component_templates::find(previous_component.name);
        const auto* const current_template = component_templates::find(current_component.name);
        if (previous_template == nullptr || current_template == nullptr ||
            previous_template->kind != kind || current_template->kind != kind ||
            previous.chain_name != current.chain_name ||
            !sequentially_adjacent(*previous.residue, *current.residue)) {
            continue;
        }

        const auto first = previous_component.find_atom(previous_atom);
        const auto second = current_component.find_atom(current_atom);
        if (first.has_value() && second.has_value()) {
            bonds.add(*first, *second, core::BondOrder::SINGLE);
        }
    }
}

[[nodiscard]] auto bond_order(const std::string& value) -> std::optional<core::BondOrder> {
    try {
        switch (::gemmi::bond_type_from_string(value)) {
        case ::gemmi::BondType::Single:
            return core::BondOrder::SINGLE;
        case ::gemmi::BondType::Double:
            return core::BondOrder::DOUBLE;
        case ::gemmi::BondType::Triple:
            return core::BondOrder::TRIPLE;
        case ::gemmi::BondType::Aromatic:
            return core::BondOrder::SINGLE;
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

    for (const auto& residue : residues) {
        const auto component = residue.component();
        const auto* const component_template = component_templates::find(component.name);
        if (component_template == nullptr) {
            continue;
        }

        add_component_bonds(result, component, component_template->bonds);
        if (component_template->kind == component_templates::ComponentKind::amino_acid) {
            if (const auto nitrogen = component.find_atom("N")) {
                // Supplement the base template with the remaining terminal hydrogen names.
                for (const auto name : {component.name == "PRO" ? "H2" : "H1", "H3"}) {
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
    -> std::vector<core::Bond> {
    if (structure.models.empty()) {
        return {};
    }

    BondAccumulator result;
    add_structure_connections(result, structure, model);

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

    return result.take();
}

auto explicit_mmcif(const ::gemmi::Structure& structure, ::gemmi::cif::Block& block,
                    const selection::SelectedModel& model) -> std::vector<core::Bond> {
    if (structure.models.empty()) {
        return {};
    }

    const auto& residues = model.residues();
    BondAccumulator result;

    struct ParsedComponentBond {
        std::string first;
        std::string second;
        core::BondOrder order;
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
    add_structure_connections(result, structure, model, connection_orders);

    return result.take();
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
