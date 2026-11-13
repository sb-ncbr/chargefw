#include "io/native/mol2_atom_types.h"

#include <chargefw/core/bond.h>
#include <chargefw/core/periodic_table.h>

#include <algorithm>
#include <cstddef>
#include <ranges>
#include <string_view>
#include <vector>

namespace chargefw::io::native::detail {
namespace {

constexpr auto hydrogen = 1;
constexpr auto carbon = 6;
constexpr auto nitrogen = 7;
constexpr auto oxygen = 8;
constexpr auto phosphorus = 15;
constexpr auto sulfur = 16;

struct Neighbor {
    std::size_t atom_index;
    core::BondOrder order;
};

class Typing {
  public:
    explicit Typing(const core::Molecule& molecule)
        : molecule_{&molecule}, neighbors_(molecule.atom_count()) {
        for (const auto& bond : molecule.bonds()) {
            neighbors_[bond.first_atom_index()].push_back({bond.second_atom_index(), bond.order()});
            neighbors_[bond.second_atom_index()].push_back({bond.first_atom_index(), bond.order()});
        }
    }

    [[nodiscard]] auto type(const std::size_t atom) const -> std::string_view {
        switch (element(atom)) {
        case carbon:
            return carbon_type(atom);
        case nitrogen:
            return nitrogen_type(atom);
        case oxygen:
            return oxygen_type(atom);
        case sulfur:
            return sulfur_type(atom);
        case phosphorus:
            return "P.3";
        case hydrogen:
            return "H";
        default:
            return core::element_symbol(element(atom));
        }
    }

  private:
    [[nodiscard]] auto element(const std::size_t atom) const -> int {
        return molecule_->atom(atom).atomic_number();
    }

    [[nodiscard]] auto formal_charge(const std::size_t atom) const -> int {
        return molecule_->atom(atom).formal_charge();
    }

    [[nodiscard]] auto bonds_of_order(const std::size_t atom, const core::BondOrder order) const
        -> std::size_t {
        return static_cast<std::size_t>(std::ranges::count(
            neighbors_[atom], order, [](const Neighbor& neighbor) { return neighbor.order; }));
    }

    [[nodiscard]] auto has_multiple_bond(const std::size_t atom) const -> bool {
        return bonds_of_order(atom, core::BondOrder::SINGLE) != neighbors_[atom].size();
    }

    [[nodiscard]] auto is_terminal_oxygen(const std::size_t atom) const -> bool {
        return element(atom) == oxygen && neighbors_[atom].size() == 1;
    }

    // A carbon with a double bond to oxygen or sulfur, as in amides and thioamides.
    [[nodiscard]] auto is_carbonyl_carbon(const std::size_t atom) const -> bool {
        return element(atom) == carbon &&
               std::ranges::any_of(neighbors_[atom], [this](const Neighbor& neighbor) {
                   return neighbor.order == core::BondOrder::DOUBLE &&
                          (element(neighbor.atom_index) == oxygen ||
                           element(neighbor.atom_index) == sulfur);
               });
    }

    // A carboxylate carbon or a phosphate phosphorus, the groups whose terminal oxygens MOL2 types
    // as O.co2: one terminal oxygen double-bonded and neutral, the others single-bonded with charge
    // -1, two terminal oxygens on carbon. Readers such as RDKit rebuild exactly this pattern from
    // the types, so carbonate and protonated groups keep O.2 and O.3.
    [[nodiscard]] auto is_carboxylate_or_phosphate(const std::size_t atom) const -> bool {
        if (element(atom) != carbon && element(atom) != phosphorus) {
            return false;
        }
        auto oxo = std::size_t{0};
        auto oxide = std::size_t{0};
        for (const auto& neighbor : neighbors_[atom]) {
            const auto partner = neighbor.atom_index;
            if (!is_terminal_oxygen(partner)) {
                continue;
            }
            if (neighbor.order == core::BondOrder::DOUBLE && formal_charge(partner) == 0) {
                ++oxo;
            } else if (neighbor.order == core::BondOrder::SINGLE && formal_charge(partner) == -1) {
                ++oxide;
            } else {
                return false;
            }
        }
        return oxo == 1 && oxide >= 1 && (element(atom) == phosphorus || oxide == 1);
    }

    // An oxygen bound to sulfur only, by a double bond or, with a positive sulfur, as an oxide:
    // S=O or S+-O-.
    [[nodiscard]] auto is_sulfur_oxo(const std::size_t sulfur_atom, const Neighbor& neighbor) const
        -> bool {
        const auto partner = neighbor.atom_index;
        if (!is_terminal_oxygen(partner)) {
            return false;
        }
        return neighbor.order == core::BondOrder::DOUBLE ||
               (formal_charge(partner) < 0 && formal_charge(sulfur_atom) > 0);
    }

    [[nodiscard]] auto carbon_type(const std::size_t atom) const -> std::string_view {
        const auto doubles = bonds_of_order(atom, core::BondOrder::DOUBLE);
        if (bonds_of_order(atom, core::BondOrder::TRIPLE) > 0 || doubles >= 2) {
            return "C.1";
        }
        if (doubles == 0) {
            return "C.3";
        }
        // Amidinium and guanidinium: C=N+ with another nitrogen.
        const auto& bonded = neighbors_[atom];
        const auto iminium = std::ranges::any_of(bonded, [this](const Neighbor& neighbor) {
            return neighbor.order == core::BondOrder::DOUBLE &&
                   element(neighbor.atom_index) == nitrogen &&
                   formal_charge(neighbor.atom_index) > 0;
        });
        const auto nitrogens = std::ranges::count_if(bonded, [this](const Neighbor& neighbor) {
            return element(neighbor.atom_index) == nitrogen;
        });
        return iminium && nitrogens >= 2 ? "C.cat" : "C.2";
    }

    [[nodiscard]] auto nitrogen_type(const std::size_t atom) const -> std::string_view {
        const auto& bonded = neighbors_[atom];
        const auto doubles = bonds_of_order(atom, core::BondOrder::DOUBLE);
        if (bonds_of_order(atom, core::BondOrder::TRIPLE) > 0) {
            return "N.1";
        }
        // A charged or pentavalent nitrogen with a double bond is trigonal planar: iminium,
        // pyridinium, and nitro nitrogens.
        if (doubles >= 2) {
            return bonded.size() >= 3 ? "N.pl3" : "N.1";
        }
        if (doubles == 1) {
            return formal_charge(atom) > 0 ? "N.pl3" : "N.2";
        }
        if (bonded.size() == 4 && formal_charge(atom) > 0) {
            return "N.4";
        }
        if (std::ranges::any_of(bonded, [this](const Neighbor& neighbor) {
                return is_carbonyl_carbon(neighbor.atom_index);
            })) {
            return "N.am";
        }
        // Conjugated with a multiple bond of carbon or nitrogen, as in anilines, enamines, and
        // guanidines.
        if (std::ranges::any_of(bonded, [this](const Neighbor& neighbor) {
                const auto partner = neighbor.atom_index;
                return (element(partner) == carbon || element(partner) == nitrogen) &&
                       has_multiple_bond(partner);
            })) {
            return "N.pl3";
        }
        return "N.3";
    }

    [[nodiscard]] auto oxygen_type(const std::size_t atom) const -> std::string_view {
        const auto& bonded = neighbors_[atom];
        if (bonded.size() == 1) {
            const auto partner = bonded.front().atom_index;
            if (is_carboxylate_or_phosphate(partner)) {
                return "O.co2";
            }
        }
        return bonds_of_order(atom, core::BondOrder::DOUBLE) > 0 ? "O.2" : "O.3";
    }

    [[nodiscard]] auto sulfur_type(const std::size_t atom) const -> std::string_view {
        const auto& bonded = neighbors_[atom];
        // Sulfoxides and sulfones in either notation; a sulfenate, S-O-, is not one.
        const auto oxo = std::ranges::count_if(bonded, [this, atom](const Neighbor& neighbor) {
            return is_sulfur_oxo(atom, neighbor);
        });
        if (oxo >= 2) {
            return "S.O2";
        }
        if (oxo == 1) {
            return "S.O";
        }
        return bonds_of_order(atom, core::BondOrder::DOUBLE) > 0 ? "S.2" : "S.3";
    }

    const core::Molecule* molecule_;
    std::vector<std::vector<Neighbor>> neighbors_;
};

} // namespace

auto sybyl_atom_types(const core::Molecule& molecule) -> std::vector<std::string_view> {
    const auto typing = Typing{molecule};
    auto types = std::vector<std::string_view>{};
    types.reserve(molecule.atom_count());
    for (std::size_t atom = 0; atom < molecule.atom_count(); ++atom) {
        types.push_back(typing.type(atom));
    }
    return types;
}

} // namespace chargefw::io::native::detail
