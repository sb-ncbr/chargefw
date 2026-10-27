#include <chargefw/charges/charge_collection.h>
#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/molecule_collection.h>

#include "calculation/fixed_charge_partition.h"

#include <snitch/snitch.hpp>

#include <algorithm>
#include <ranges>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace calculation = chargefw::calculation;
namespace core = chargefw::core;

namespace {

struct PartitionInput {
    core::MoleculeCollection molecules;
    calculation::FixedIons fixed_ions;
};

auto make_partition_input(const std::vector<calculation::FixedAtomCharge>& sources)
    -> PartitionInput {
    auto first_affected = core::Molecule{
        std::vector{core::Atom{6, 1, "C0"}, core::Atom{7, 0, "N1"}, core::Atom{12, 2, "Mg2"},
                    core::Atom{8, -1, "O3"}},
        std::vector{core::Bond{3, 0, core::BondOrder::TRIPLE},
                    core::Bond{0, 1, core::BondOrder::SINGLE}},
        {core::Conformer{{{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {3.0, 0.0, 0.0}, {2.0, 1.0, 0.0}},
                         "alpha"},
         core::Conformer{{{0.0, 0.0, 1.0}, {1.0, 0.0, 1.0}, {3.0, 0.0, 1.0}, {2.0, 1.0, 1.0}},
                         "beta"}},
        "affected-one"};
    auto second_affected = core::Molecule{
        std::vector{core::Atom{1, 1, "H0"}, core::Atom{6, 0, "C1"}, core::Atom{7, 0, "N2"},
                    core::Atom{8, -1, "O3"}},
        std::vector{core::Bond{2, 1, core::BondOrder::DOUBLE}},
        {core::Conformer{{{5.0, 0.0, 0.0}, {6.0, 0.0, 0.0}, {7.0, 0.0, 0.0}, {8.0, 0.0, 0.0}},
                         "gamma"},
         core::Conformer{{{5.0, 0.0, 1.0}, {6.0, 0.0, 1.0}, {7.0, 0.0, 1.0}, {8.0, 0.0, 1.0}},
                         "delta"}},
        "affected-three"};
    auto molecules = core::MoleculeCollection{
        std::vector{
            core::Molecule{std::vector{core::Atom{6, 1, "before-C"}, core::Atom{8, -1, "before-O"}},
                           std::vector{core::Bond{1, 0, core::BondOrder::DOUBLE}},
                           {},
                           "unaffected-before"},
            std::move(first_affected),
            core::Molecule{
                std::vector{core::Atom{9, 0, "between-F"}}, {}, {}, "unaffected-between"},
            std::move(second_affected)},
        "named-collection"};
    return {std::move(molecules), calculation::FixedIons{sources}};
}

auto make_partition_from_temporary_inputs() -> calculation::detail::FixedChargePartition {
    auto input = make_partition_input({{3, 3, -0.25}, {1, 2, 0.5}, {3, 0, 0.75}});
    return calculation::detail::make_fixed_charge_partition(input.molecules, input.fixed_ions);
}

auto make_active_charge_set() -> chargefw::charges::ChargeSet {
    using chargefw::charges::AtomicCharges;
    using chargefw::charges::ChargeAssignment;
    using chargefw::charges::ChargeTarget;
    return chargefw::charges::ChargeSet{
        "test-method",
        std::vector<ChargeAssignment>{
            {ChargeTarget{3, 1}, AtomicCharges{{10.0, 20.0}}},
            {ChargeTarget{0, std::nullopt}, AtomicCharges{{7.0, 8.0}}},
            {ChargeTarget{1, 0}, AtomicCharges{{-1.0, -2.0, -3.0}}},
            {ChargeTarget{2, std::nullopt}, AtomicCharges{{-4.0}}},
            {ChargeTarget{3, 0}, AtomicCharges{{-10.0, -20.0}}},
            {ChargeTarget{1, 1}, AtomicCharges{{1.0, 2.0, 3.0}}},
        },
        "test-parameters"};
}

} // namespace

TEST_CASE("fixed-charge partition preserves original ordering and owned mappings",
          "[calculation][fixed-charge-partition]") {
    const auto input = make_partition_input({{3, 3, -0.25}, {1, 2, 0.5}, {3, 0, 0.75}});
    const auto partition =
        calculation::detail::make_fixed_charge_partition(input.molecules, input.fixed_ions);

    struct ExpectedTarget {
        std::vector<std::size_t> active_atoms;
        std::vector<std::size_t> active_bonds;
        // Sources sorted by original atom index: (atom index, charge).
        std::vector<std::pair<std::size_t, double>> sources;
    };
    const auto expected = std::vector<ExpectedTarget>{
        {.active_atoms = {0, 1}, .active_bonds = {0}},
        {.active_atoms = {0, 1, 3}, .active_bonds = {0, 1}, .sources = {{2, 0.5}}},
        {.active_atoms = {0}},
        {.active_atoms = {1, 2}, .active_bonds = {0}, .sources = {{0, 0.75}, {3, -0.25}}}};

    CHECK(partition.active_molecules.name() == "named-collection");
    REQUIRE(partition.active_molecules.size() == expected.size());
    REQUIRE(partition.targets.size() == expected.size());
    const auto same_position = [](const core::Position& a, const core::Position& b) {
        return a.x == b.x && a.y == b.y && a.z == b.z;
    };

    for (std::size_t index = 0; index < expected.size(); ++index) {
        CAPTURE(index);
        const auto& target = partition.targets[index];
        const auto& original = input.molecules[index];
        const auto& active = partition.active_molecules[index];
        CHECK(target.active_atom_indices == expected[index].active_atoms);
        CHECK(target.active_bond_indices == expected[index].active_bonds);
        CHECK(target.active_charge == 0.0);
        CHECK(core::total_formal_charge(active) == target.active_charge);

        // Fixed sources keep their original identity and conformer-local positions.
        REQUIRE(target.sources.size() == expected[index].sources.size());
        for (std::size_t source = 0; source < target.sources.size(); ++source) {
            CHECK(target.sources[source].molecule_index == index);
            CHECK(target.sources[source].atom_index == expected[index].sources[source].first);
            CHECK(target.sources[source].charge == expected[index].sources[source].second);
        }
        REQUIRE(target.source_positions.size() ==
                (target.sources.empty() ? 0 : original.conformer_count()));
        for (std::size_t conformer = 0; conformer < target.source_positions.size(); ++conformer) {
            REQUIRE(target.source_positions[conformer].size() == target.sources.size());
            for (std::size_t source = 0; source < target.sources.size(); ++source) {
                CHECK(same_position(
                    target.source_positions[conformer][source],
                    original.conformer(conformer).positions()[target.sources[source].atom_index]));
            }
        }

        // Active molecules are the original atoms, bonds, and conformers in active order.
        CHECK(active.name() == original.name());
        REQUIRE(active.atom_count() == target.active_atom_indices.size());
        for (std::size_t atom = 0; atom < active.atom_count(); ++atom) {
            const auto& source_atom = original.atom(target.active_atom_indices[atom]);
            CHECK(active.atom(atom).name() == source_atom.name());
            CHECK(active.atom(atom).formal_charge() == source_atom.formal_charge());
        }
        REQUIRE(active.bond_count() == target.active_bond_indices.size());
        for (std::size_t bond = 0; bond < active.bond_count(); ++bond) {
            const auto& source_bond = original.bond(target.active_bond_indices[bond]);
            CHECK(target.active_atom_indices[active.bond(bond).first_atom_index()] ==
                  source_bond.first_atom_index());
            CHECK(target.active_atom_indices[active.bond(bond).second_atom_index()] ==
                  source_bond.second_atom_index());
            CHECK(active.bond(bond).order() == source_bond.order());
        }
        REQUIRE(active.conformer_count() == original.conformer_count());
        for (std::size_t conformer = 0; conformer < active.conformer_count(); ++conformer) {
            CHECK(active.conformer(conformer).name() == original.conformer(conformer).name());
            for (std::size_t atom = 0; atom < active.atom_count(); ++atom) {
                CHECK(same_position(
                    active.conformer(conformer).positions()[atom],
                    original.conformer(conformer).positions()[target.active_atom_indices[atom]]));
            }
        }
    }
}

TEST_CASE("fixed-charge partition owns data after input destruction and moves",
          "[calculation][fixed-charge-partition]") {
    auto partition = make_partition_from_temporary_inputs();
    auto moved_partition = std::move(partition);

    CHECK(moved_partition.active_molecules.name() == "named-collection");
    CHECK(moved_partition.active_molecules[1].atom(2).name() == "O3");
    CHECK(moved_partition.targets[1].sources[0].charge == 0.5);
    CHECK(moved_partition.targets[3].source_positions[1][1].z == 1.0);
}

TEST_CASE("fixed-charge partition retains disconnected active components as one target",
          "[calculation][fixed-charge-partition]") {
    const auto molecules = core::MoleculeCollection{std::vector{core::Molecule{
        std::vector{core::Atom{6, 1}, core::Atom{12, 2}, core::Atom{6, 3}, core::Atom{8, 4}},
        std::vector{core::Bond{3, 2}},
        {core::Conformer{{{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {2.0, 0.0, 0.0}, {3.0, 0.0, 0.0}}}},
        "disconnected-active"}}};
    const auto fixed_ions = calculation::FixedIons{{{0, 1, 0.75}}};
    const auto partition = calculation::detail::make_fixed_charge_partition(molecules, fixed_ions);

    REQUIRE(partition.active_molecules.size() == 1);
    REQUIRE(partition.targets.size() == 1);
    CHECK(partition.targets[0].active_atom_indices == std::vector<std::size_t>{0, 2, 3});
    CHECK(partition.targets[0].active_charge == 8.0);
    CHECK(core::total_formal_charge(partition.active_molecules[0]) == 8.0);
    CHECK(partition.active_molecules[0].atom_count() == 3);
    REQUIRE(partition.active_molecules[0].bond_count() == 1);
    CHECK(partition.active_molecules[0].bond(0).first_atom_index() == 2);
    CHECK(partition.active_molecules[0].bond(0).second_atom_index() == 1);
}

TEST_CASE("empty fixed-charge partition is an identity copy",
          "[calculation][fixed-charge-partition]") {
    auto input = make_partition_input({});
    const auto partition =
        calculation::detail::make_fixed_charge_partition(input.molecules, input.fixed_ions);

    REQUIRE(partition.active_molecules.size() == input.molecules.size());
    REQUIRE(partition.targets.size() == input.molecules.size());
    for (std::size_t index = 0; index < input.molecules.size(); ++index) {
        CHECK(partition.targets[index].sources.empty());
        CHECK(partition.targets[index].source_positions.empty());
        CHECK(partition.targets[index].active_charge ==
              core::total_formal_charge(input.molecules[index]));
        CHECK(partition.targets[index].active_atom_indices.size() ==
              input.molecules[index].atom_count());
        CHECK(partition.targets[index].active_bond_indices.size() ==
              input.molecules[index].bond_count());
        for (std::size_t atom_index = 0; atom_index < input.molecules[index].atom_count();
             ++atom_index) {
            CHECK(partition.targets[index].active_atom_indices[atom_index] == atom_index);
        }
        for (std::size_t bond_index = 0; bond_index < input.molecules[index].bond_count();
             ++bond_index) {
            CHECK(partition.targets[index].active_bond_indices[bond_index] == bond_index);
        }
    }
}

TEST_CASE("fixed-charge reassembly scatters active charges and preserves assignment metadata",
          "[calculation][fixed-charge-partition]") {
    const auto result = [] {
        auto input = make_partition_input({{3, 3, -0.25}, {1, 2, 0.5}, {3, 0, 0.75}});
        const auto partition =
            calculation::detail::make_fixed_charge_partition(input.molecules, input.fixed_ions);
        const auto active_charges = make_active_charge_set();
        const auto reassembled =
            calculation::detail::reassemble_fixed_charge_results(active_charges, partition);

        return reassembled;
    }();

    CHECK(result.method_id() == "test-method");
    REQUIRE(result.parameter_set_id().has_value());
    CHECK(*result.parameter_set_id() == "test-parameters");
    REQUIRE(result.size() == 6);

    const auto expected = chargefw::charges::ChargeSet{
        "test-method",
        {{chargefw::charges::ChargeTarget{3, 1},
          chargefw::charges::AtomicCharges{{0.75, 10.0, 20.0, -0.25}}},
         {chargefw::charges::ChargeTarget{0, std::nullopt},
          chargefw::charges::AtomicCharges{{7.0, 8.0}}},
         {chargefw::charges::ChargeTarget{1, 0},
          chargefw::charges::AtomicCharges{{-1.0, -2.0, 0.5, -3.0}}},
         {chargefw::charges::ChargeTarget{2, std::nullopt},
          chargefw::charges::AtomicCharges{{-4.0}}},
         {chargefw::charges::ChargeTarget{3, 0},
          chargefw::charges::AtomicCharges{{0.75, -10.0, -20.0, -0.25}}},
         {chargefw::charges::ChargeTarget{1, 1},
          chargefw::charges::AtomicCharges{{1.0, 2.0, 0.5, 3.0}}}},
    };
    for (std::size_t index = 0; index < expected.size(); ++index) {
        CAPTURE(index);
        const auto& actual = result.assignment(index);
        const auto& reference = expected.assignment(index);
        CHECK(actual.target.molecule_index == reference.target.molecule_index);
        CHECK(actual.target.conformer_index == reference.target.conformer_index);
        CHECK(std::ranges::equal(actual.charges.values(), reference.charges.values()));
    }
    CHECK(result.assignment(0).charges.total() == 30.5);
    CHECK(result.assignment(2).charges.total() == -5.5);
}

TEST_CASE("fixed-charge reassembly supports identity partitions and absent parameter IDs",
          "[calculation][fixed-charge-partition]") {
    auto input = make_partition_input({});
    const auto partition =
        calculation::detail::make_fixed_charge_partition(input.molecules, input.fixed_ions);
    const auto active = chargefw::charges::ChargeSet{
        "identity-method",
        {{chargefw::charges::ChargeTarget{1, std::nullopt},
          chargefw::charges::AtomicCharges{{1.0, 2.0, 3.0, 4.0}}}},
    };
    const auto result = calculation::detail::reassemble_fixed_charge_results(active, partition);

    CHECK(result.method_id() == "identity-method");
    CHECK_FALSE(result.parameter_set_id().has_value());
    REQUIRE(result.size() == 1);
    CHECK(result.assignment(0).target.molecule_index == 1);
    CHECK_FALSE(result.assignment(0).target.conformer_index.has_value());
    CHECK(std::ranges::equal(result.assignment(0).charges.values(),
                             std::vector<double>{1.0, 2.0, 3.0, 4.0}));
}

TEST_CASE("fixed-charge reassembly rejects target and active-size mismatches",
          "[calculation][fixed-charge-partition]") {
    auto input = make_partition_input({{1, 2, 0.5}});
    const auto partition =
        calculation::detail::make_fixed_charge_partition(input.molecules, input.fixed_ions);
    const auto wrong_size = chargefw::charges::ChargeSet{
        "test-method",
        {{chargefw::charges::ChargeTarget{1, std::nullopt},
          chargefw::charges::AtomicCharges{{1.0, 2.0}}}},
    };
    CHECK_THROWS_AS(calculation::detail::reassemble_fixed_charge_results(wrong_size, partition),
                    std::invalid_argument);

    const auto wrong_target = chargefw::charges::ChargeSet{
        "test-method",
        {{chargefw::charges::ChargeTarget{99, std::nullopt}, chargefw::charges::AtomicCharges{{}}}},
    };
    CHECK_THROWS_AS(calculation::detail::reassemble_fixed_charge_results(wrong_target, partition),
                    std::out_of_range);
}

TEST_CASE("fixed-charge partition factory shares source validation",
          "[calculation][fixed-charge-partition]") {
    const auto unbonded = core::MoleculeCollection{
        std::vector{core::Molecule{std::vector{core::Atom{6}, core::Atom{6}},
                                   {},
                                   {core::Conformer{{{0.0, 0.0, 0.0}, {2.0, 0.0, 0.0}}}}}}};
    const auto bonded = core::MoleculeCollection{
        std::vector{core::Molecule{std::vector{core::Atom{6}, core::Atom{6}},
                                   std::vector{core::Bond{0, 1}},
                                   {core::Conformer{{{0.0, 0.0, 0.0}, {2.0, 0.0, 0.0}}}}}}};
    const auto fixed_ions = calculation::FixedIons{{{0, 0, 1.0}}};

    CHECK_THROWS_AS(calculation::detail::make_fixed_charge_partition(
                        unbonded, calculation::FixedIons{{{1, 0, 1.0}}}),
                    std::invalid_argument);
    CHECK_THROWS_AS(calculation::detail::make_fixed_charge_partition(bonded, fixed_ions),
                    std::invalid_argument);
    const auto coincident = core::MoleculeCollection{
        std::vector{core::Molecule{std::vector{core::Atom{6}, core::Atom{6}},
                                   {},
                                   {core::Conformer{{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}}}}}};
    CHECK_THROWS_AS(calculation::detail::make_fixed_charge_partition(coincident, fixed_ions),
                    std::invalid_argument);
}
