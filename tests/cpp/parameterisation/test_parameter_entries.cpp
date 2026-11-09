#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <chargefw/calculation/assessment.h>
#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/core/position.h>
#include <chargefw/parameterisation/parameter_entries.h>
#include <chargefw/parameters/models/atom_parameters.h>
#include <chargefw/parameters/models/bond_parameters.h>
#include <chargefw/parameters/models/parameter_set.h>
#include <chargefw/parameters/models/parameter_set_metadata.h>

#include <cstddef>
#include <snitch/snitch.hpp>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// What the assessed molecules show about each parameter entry. The SQE sets have a plain hydrogen
// entry 0 and oxygen entry 1.

namespace calculation = chargefw::calculation;
namespace core = chargefw::core;
namespace parameterisation = chargefw::parameterisation;
namespace parameters = chargefw::parameters;

namespace {

auto molecule(std::vector<core::Atom> atoms, std::vector<core::Bond> bonds,
              std::vector<core::Position> positions, std::string name) -> core::Molecule {
    auto conformers = std::vector{core::Conformer{std::move(positions), "model-1"}};
    return core::Molecule{std::move(atoms), std::move(bonds), std::move(conformers),
                          std::move(name)};
}

auto hydrogen_peroxide() -> core::Molecule {
    return molecule({core::Atom{1, 0}, core::Atom{8, 0}, core::Atom{8, 0}, core::Atom{1, 0}},
                    {core::Bond{0, 1}, core::Bond{1, 2}, core::Bond{2, 3}},
                    {{.x = -0.9, .y = 0.3, .z = 0.0},
                     {.x = 0.0, .y = 0.0, .z = 0.0},
                     {.x = 1.47, .y = 0.0, .z = 0.0},
                     {.x = 2.37, .y = 0.3, .z = 0.4}},
                    "hydrogen-peroxide");
}

// O-H: both atoms have one bond.
auto hydroxyl() -> core::Molecule {
    return molecule({core::Atom{8, 0}, core::Atom{1, 0}}, {core::Bond{0, 1}},
                    {{.x = 0.0, .y = 0.0, .z = 0.0}, {.x = 0.97, .y = 0.0, .z = 0.0}}, "hydroxyl");
}

auto sqe_set(std::vector<parameters::BondParameterEntry> bonds) -> parameters::ParameterSet {
    const auto values = [](const double electronegativity) {
        return std::vector<parameters::NamedParameter>{
            {.name = "electronegativity", .value = electronegativity},
            {.name = "hardness", .value = 13.0},
            {.name = "width", .value = 1.0}};
    };
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{.id = "entries", .method_id = "sqe"},
        {},
        parameters::AtomParameters{
            {{.key = chargefw::test::plain_atom_key(1), .parameters = values(4.5)},
             {.key = chargefw::test::plain_atom_key(8), .parameters = values(8.7)}}},
        parameters::BondParameters{std::move(bonds)}};
}

auto kappa(parameters::BondParameterKey key) -> parameters::BondParameterEntry {
    return {.key = std::move(key), .parameters = {{.name = "kappa", .value = 1.0}}};
}

auto assess(std::vector<core::Molecule> molecules, parameters::ParameterSet set)
    -> calculation::AssessmentResult {
    auto result = calculation::assess(
        calculation::AssessmentRequest{.molecules = core::MoleculeCollection{std::move(molecules)},
                                       .parameter_sets = {std::move(set)},
                                       .method_id = "sqe",
                                       .parameter_set_id = "entries"});
    REQUIRE(result.default_plan() != nullptr);
    return result;
}

} // namespace

TEST_CASE("Entries report their atoms, bonds, and the bonds that follow single-bonded atoms",
          "[parameterisation][parameter-entries]") {
    const auto molecules = std::vector{chargefw::test::make_water(), hydrogen_peroxide()};

    // H-O and O-O entries: every hydrogen has one bond, and only hydrogens end the H-O bonds.
    const auto separate = assess(molecules, sqe_set({kappa(chargefw::test::plain_bond_key(1, 8)),
                                                     kappa(chargefw::test::plain_bond_key(8, 8))}));
    const auto counts = parameterisation::entry_counts(separate, *separate.default_plan());
    CHECK(counts.atoms == std::vector<std::size_t>{4, 3});
    CHECK(counts.bonds == std::vector<std::size_t>{4, 1});
    const auto entries = parameterisation::single_bond_entries(separate, *separate.default_plan());
    REQUIRE(entries.size() == 2);
    CHECK(entries[0].single_bonded);
    CHECK(entries[0].bond_entries == std::vector<std::size_t>{0});
    CHECK(entries[0].exclusive_bonds);
    CHECK_FALSE(entries[1].single_bonded);

    // One entry for every bond: the O-O bond has no single-bonded atom, so the entry does not
    // follow the hydrogens alone.
    const auto shared = assess(molecules, sqe_set({kappa(chargefw::test::plain_bond_key())}));
    const auto shared_entries =
        parameterisation::single_bond_entries(shared, *shared.default_plan());
    CHECK(shared_entries[0].single_bonded);
    CHECK(shared_entries[0].bond_entries == std::vector<std::size_t>{0});
    CHECK_FALSE(shared_entries[0].exclusive_bonds);

    // In O-H both atoms have one bond, so neither entry is single-bonded, though water's hydrogens
    // alone would be.
    const auto radical = assess({chargefw::test::make_water(), hydroxyl()},
                                sqe_set({kappa(chargefw::test::plain_bond_key(1, 8))}));
    for (const auto& entry :
         parameterisation::single_bond_entries(radical, *radical.default_plan())) {
        CHECK_FALSE(entry.single_bonded);
    }
    // Also for a method without bond parameters, whose classification lists no bonds.
    const auto eem_radical = calculation::assess(calculation::AssessmentRequest{
        .molecules =
            core::MoleculeCollection{std::vector{chargefw::test::make_water(), hydroxyl()}},
        .parameter_sets = {chargefw::test::make_eem_ho_parameters(0.5)},
        .method_id = "eem",
        .parameter_set_id = "test-eem-ho"});
    REQUIRE(eem_radical.default_plan() != nullptr);
    for (const auto& entry :
         parameterisation::single_bond_entries(eem_radical, *eem_radical.default_plan())) {
        CHECK_FALSE(entry.single_bonded);
        CHECK(entry.bond_entries.empty());
    }
}

TEST_CASE("Entries count the calculated atoms, not fixed ions",
          "[parameterisation][parameter-entries]") {
    // Water (O, H, H) and a magnesium ion with a fixed charge; the H/O set has no entry for it.
    auto atoms =
        std::vector{core::Atom{8, 0}, core::Atom{1, 0}, core::Atom{1, 0}, core::Atom{12, 2}};
    auto positions = std::vector<core::Position>{
        {0.0, 0.0, 0.0}, {0.96, 0.0, 0.0}, {-0.24, 0.93, 0.0}, {5.0, 0.0, 0.0}};
    const auto water_with_ion = molecule(std::move(atoms), {core::Bond{0, 1}, core::Bond{0, 2}},
                                         std::move(positions), "water-with-ion");
    const auto assessment = calculation::assess(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{water_with_ion}},
        .parameter_sets = {chargefw::test::make_eem_ho_parameters(0.5)},
        .method_id = "eem",
        .parameter_set_id = "test-eem-ho",
        .fixed_ions = calculation::FixedIons{
            .sources = {{.molecule_index = 0, .atom_index = 3, .charge = 2.0}}}});
    const auto* plan = assessment.default_plan();
    REQUIRE(plan != nullptr);

    const auto counts = parameterisation::entry_counts(assessment, *plan);
    CHECK(counts.atoms == std::vector<std::size_t>{2, 1});
    CHECK(counts.bonds.empty());
    // The calculated part does not line up with the molecule, so no entry is single-bonded.
    for (const auto& entry : parameterisation::single_bond_entries(assessment, *plan)) {
        CHECK_FALSE(entry.single_bonded);
    }

    const auto formal = calculation::assess(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .method_id = "formal"});
    REQUIRE(formal.default_plan() != nullptr);
    CHECK_THROWS_AS(parameterisation::entry_counts(formal, *formal.default_plan()),
                    std::invalid_argument);
    CHECK_THROWS_AS(parameterisation::single_bond_entries(assessment, *formal.default_plan()),
                    std::invalid_argument);
}
