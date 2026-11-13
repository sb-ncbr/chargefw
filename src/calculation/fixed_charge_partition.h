#pragma once

#include <chargefw/calculation/assessment.h>
#include <chargefw/charges/charge_set.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/calculation_input.h>

#include <cstddef>
#include <optional>
#include <vector>

namespace chargefw::calculation::detail {

struct FixedChargePartitionTarget {
    // Active indices map to atoms/bonds in the original molecule at this target index.
    std::vector<std::size_t> active_atom_indices;
    std::vector<std::size_t> active_bond_indices;
    std::vector<FixedAtomCharge> sources;
    // Affected targets use original conformer order and sorted source atom order.
    std::vector<std::vector<core::Position>> source_positions;
    double active_charge = 0.0;
};

struct FixedChargePartition {
    core::MoleculeCollection active_molecules;
    std::vector<FixedChargePartitionTarget> targets;
};

[[nodiscard]] auto make_fixed_charge_partition(const core::MoleculeCollection& molecules,
                                               const FixedIons& fixed_ions) -> FixedChargePartition;

auto validate_partition_active_molecules(const features::PreparedMoleculeCollection& molecules,
                                         const FixedChargePartition& partition) -> void;

[[nodiscard]] auto materialize_fixed_point_sources(const FixedChargePartitionTarget& target,
                                                   std::optional<std::size_t> conformer_index)
    -> std::vector<methods::FixedPointSource>;

// The partition must come from the validated factory and active assignments must match its active
// molecules.
[[nodiscard]] auto reassemble_fixed_charge_results(const charges::ChargeSet& active_charges,
                                                   const FixedChargePartition& partition)
    -> charges::ChargeSet;

} // namespace chargefw::calculation::detail
