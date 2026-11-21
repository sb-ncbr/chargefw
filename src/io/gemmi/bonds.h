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

// On input, unknown_order_indices lists ascending indices into explicit_bonds (PDB only).
// On return, it lists the corresponding bonds still lacking orders, indexed in the result.
// residues_with_file_bonds is ExplicitBonds::residues_with_file_bonds, or empty when the file's
// bonds were not read.
// The dictionary, when given, adds the bonds within residues that have no built-in template; under
// hybrid, only within those the file states no bond within either. Under hybrid it also gives the
// orders of the pairs listed in unknown_order_indices, in a residue without a built-in template
// whose entry has the same bonds between non-hydrogen atoms as the file. An entry that contradicts
// an element of the residue is not used at all, and a residue with a temporary ligand name, such as
// LIG, is not looked up. dictionary_issue receives the warning that names the residues that took
// nothing from the dictionary for one of these reasons.
[[nodiscard]] auto assign(const selection::SelectedModel& model, BondStrategy strategy,
                          const ComponentDictionary* dictionary,
                          std::vector<core::Bond> explicit_bonds,
                          const std::vector<bool>& residues_with_file_bonds,
                          ExplicitBondOrders explicit_orders,
                          std::vector<std::size_t>& unknown_order_indices,
                          std::optional<ImportIssue>& dictionary_issue) -> std::vector<core::Bond>;

// Bonds stated in the file. Connections to atoms of a symmetry-related copy are not bonds of
// the molecule; they are counted, not added.
struct ExplicitBonds {
    std::vector<core::Bond> bonds;
    std::size_t skipped_symmetry_connections = 0;
    // Indices in bonds of the aromatic bonds, read as single bonds. The record keeps explicit bonds
    // first and in this order, so these are its bond indices too.
    std::vector<std::size_t> aromatic_bonds = {};
    // Intra-residue pairs of non-hydrogen atoms stated in PDB CONECT. Assignment removes
    // template-covered pairs and
    // maps these explicit bond indices to the retained topology.
    std::vector<std::size_t> unknown_order_indices = {};
    // For each of the model's selected residues, whether the file states a bond within it, read or
    // not: a _chem_comp_bond row counts for every residue of its component, also when it names an
    // absent atom or an order that is not read. A connection to another residue does not count.
    std::vector<bool> residues_with_file_bonds = {};
};

[[nodiscard]] auto explicit_pdb(const ::gemmi::Structure& structure,
                                const selection::SelectedModel& model) -> ExplicitBonds;

[[nodiscard]] auto explicit_mmcif(const ::gemmi::Structure& structure, ::gemmi::cif::Block& block,
                                  const selection::SelectedModel& model, BondStrategy strategy)
    -> ExplicitBonds;

// The record warning for skipped connections to symmetry copies.
[[nodiscard]] auto symmetry_connections_issue(std::size_t count) -> ImportIssue;

} // namespace chargefw::io::gemmi::bonds
