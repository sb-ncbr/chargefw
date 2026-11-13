#pragma once

#include <chargefw/adapters/conformer_selection.h>

#include <cstdint>
#include <string_view>

namespace chargefw::adapters::gemmi {

enum class RecordSelection : std::uint8_t {
    all,
    polymers_and_ligands,
    polymers,
};

[[nodiscard]] auto record_selection_from_string(std::string_view value) -> RecordSelection;
[[nodiscard]] auto to_string(RecordSelection value) -> std::string_view;

enum class BondStrategy : std::uint8_t {
    none,
    explicit_bonds,
    templates,
    hybrid,
};

[[nodiscard]] auto bond_strategy_from_string(std::string_view value) -> BondStrategy;
[[nodiscard]] auto to_string(BondStrategy value) -> std::string_view;

struct InputOptions {
    RecordSelection selection = RecordSelection::all;
    BondStrategy bonds = BondStrategy::none;
    ConformerSelection conformers = ConformerSelection::all;
};

} // namespace chargefw::adapters::gemmi
