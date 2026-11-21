#pragma once

#include <chargefw/io/conformer_selection.h>

#include <cstdint>
#include <memory>
#include <string_view>

namespace chargefw::io {

class ComponentDictionary;

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

struct StructureOptions {
    RecordSelection selection = RecordSelection::all;
    BondStrategy bonds = BondStrategy::none;
    ConformerSelection conformers = ConformerSelection::all;
    // Bonds within residues that the built-in templates and the file do not give, for the
    // templates and hybrid strategies. Nothing is taken from a dictionary unless one is given here.
    std::shared_ptr<const ComponentDictionary> components = {};
};

} // namespace chargefw::io
