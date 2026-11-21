#pragma once

#include <chargefw/core/bond.h>
#include <chargefw/io/molecule_record.h>
#include <chargefw/io/structure_options.h>

#include "selection.h"

#include <gemmi/model.hpp>

namespace gemmi::cif {
struct Block;
}

namespace chargefw::io {
class ComponentDictionary;
}

#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

namespace chargefw::io::gemmi::bonds {

// Whether explicit bonds carry source bond orders (mmCIF component bonds) or only connectivity
// (PDB CONECT). The hybrid strategy puts bonds with source orders before the template bonds.
enum class ExplicitBondOrders : std::uint8_t {
    connectivity_only,
    from_source,
};

// Bonds stated in the file. Connections to atoms of a symmetry-related copy are not bonds of
// the molecule; they are counted, not added.
struct ExplicitBonds {
    std::vector<core::Bond> bonds;
    ExplicitBondOrders orders = ExplicitBondOrders::from_source;
    std::size_t skipped_symmetry_connections = 0;
    // Indices in bonds of the aromatic bonds, read as single bonds. The record keeps explicit bonds
    // first and in this order, so these are its bond indices too.
    std::vector<std::size_t> aromatic_bonds = {};
    // Ascending indices in bonds of the intra-residue pairs of non-hydrogen atoms stated in PDB
    // CONECT, which gives no orders.
    std::vector<std::size_t> unknown_order_indices = {};
    // For each of the model's selected residues, whether the file states a bond within it, read or
    // not: a _chem_comp_bond row counts for every residue of its component, also when it names an
    // absent atom or an order that is not read. A connection to another residue does not count.
    std::vector<bool> residues_with_file_bonds = {};
};

struct Assignment {
    std::vector<core::Bond> bonds;
    // Ascending indices in bonds of the file's bonds that still lack an order.
    std::vector<std::size_t> unknown_order_indices = {};
    // The warning that names the residues that took nothing from the dictionary although it has an
    // entry for them, or their name is a temporary ligand name, such as LIG.
    std::optional<ImportIssue> dictionary_issue = {};
};

// The bonds of the model under a strategy. explicit_bonds are the bonds read from the file for the
// explicit and hybrid strategies.
// The dictionary, when given, adds the bonds within residues that have no built-in template; under
// hybrid, only within those the file states no bond within either. Under hybrid it also gives the
// orders of the file's bonds of unknown order, in a residue without a built-in template whose entry
// has the same bonds between non-hydrogen atoms as the file. An entry that contradicts an element
// of the residue is not used at all, and a residue with a temporary ligand name is not looked up.
[[nodiscard]] auto assign(const selection::SelectedModel& model, BondStrategy strategy,
                          const ComponentDictionary* dictionary, ExplicitBonds explicit_bonds)
    -> Assignment;

[[nodiscard]] auto explicit_pdb(const ::gemmi::Structure& structure,
                                const selection::SelectedModel& model) -> ExplicitBonds;

[[nodiscard]] auto explicit_mmcif(const ::gemmi::Structure& structure, ::gemmi::cif::Block& block,
                                  const selection::SelectedModel& model, BondStrategy strategy)
    -> ExplicitBonds;

// The record warning for skipped connections to symmetry copies.
[[nodiscard]] auto symmetry_connections_issue(std::size_t count) -> ImportIssue;

} // namespace chargefw::io::gemmi::bonds
