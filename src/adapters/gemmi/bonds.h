#pragma once

#include <chargefw/adapters/gemmi/input_options.h>
#include <chargefw/adapters/molecule_record.h>
#include <chargefw/core/bond.h>

#include "selection.h"

#include <gemmi/model.hpp>

namespace gemmi::cif {
struct Block;
}

#include <cstddef>
#include <cstdint>
#include <vector>

namespace chargefw::adapters::gemmi::bonds {

// Whether explicit bonds carry source bond orders (mmCIF component bonds) or only connectivity
// (PDB CONECT). The hybrid strategy lets source orders override template orders for the same pair.
enum class ExplicitBondOrders : std::uint8_t { connectivity_only, from_source };

[[nodiscard]] auto assign(const selection::SelectedModel& model, BondStrategy strategy,
                          std::vector<core::Bond> explicit_bonds,
                          ExplicitBondOrders explicit_orders) -> std::vector<core::Bond>;

// Bonds stated in the file. Connections to atoms of a symmetry-related copy are not bonds of
// the molecule; they are counted, not added.
struct ExplicitBonds {
    std::vector<core::Bond> bonds;
    std::size_t skipped_symmetry_connections = 0;
};

[[nodiscard]] auto explicit_pdb(const ::gemmi::Structure& structure,
                                const selection::SelectedModel& model) -> ExplicitBonds;

[[nodiscard]] auto explicit_mmcif(const ::gemmi::Structure& structure, ::gemmi::cif::Block& block,
                                  const selection::SelectedModel& model) -> ExplicitBonds;

// The record warning for skipped connections to symmetry copies.
[[nodiscard]] auto symmetry_connections_diagnostic(std::size_t count) -> MoleculeRecordDiagnostic;

} // namespace chargefw::adapters::gemmi::bonds
