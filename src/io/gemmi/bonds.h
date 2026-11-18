#pragma once

#include <chargefw/core/bond.h>
#include <chargefw/io/molecule_record.h>
#include <chargefw/io/structure_options.h>

#include "selection.h"

#include <gemmi/model.hpp>

namespace gemmi::cif {
struct Block;
}

#include <cstddef>
#include <cstdint>
#include <map>
#include <memory>
#include <span>
#include <string_view>
#include <tuple>
#include <vector>

namespace chargefw::io::gemmi::bonds {

// Whether explicit bonds carry source bond orders (mmCIF component bonds) or only connectivity
// (PDB CONECT). The hybrid strategy lets source orders override template orders for the same pair.
enum class ExplicitBondOrders : std::uint8_t {
    connectivity_only,
    from_source,
};

[[nodiscard]] auto assign(const selection::SelectedModel& model, BondStrategy strategy,
                          std::vector<core::Bond> explicit_bonds,
                          ExplicitBondOrders explicit_orders) -> std::vector<core::Bond>;

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

// The atom that a connection partner names, as the readers resolve it; null when there is none.
[[nodiscard]] auto partner_atom(const ResidueIndex& residues, const ::gemmi::AtomAddress& address)
    -> const ::gemmi::Atom*;

// Bonds stated in the file. Connections to atoms of a symmetry-related copy are not bonds of
// the molecule; they are counted, not added.
struct ExplicitBonds {
    std::vector<core::Bond> bonds;
    std::size_t skipped_symmetry_connections = 0;
    // Indices in bonds of the aromatic bonds, read as single bonds. The record keeps explicit bonds
    // first and in this order, so these are its bond indices too.
    std::vector<std::size_t> aromatic_bonds = {};
};

[[nodiscard]] auto explicit_pdb(const ::gemmi::Structure& structure,
                                const selection::SelectedModel& model) -> ExplicitBonds;

[[nodiscard]] auto explicit_mmcif(const ::gemmi::Structure& structure, ::gemmi::cif::Block& block,
                                  const selection::SelectedModel& model) -> ExplicitBonds;

// The record warning for skipped connections to symmetry copies.
[[nodiscard]] auto symmetry_connections_issue(std::size_t count) -> ImportIssue;

} // namespace chargefw::io::gemmi::bonds
