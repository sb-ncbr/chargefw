#include "support/test_charge_results.h"

#include <chargefw/charges/atomic_charges.h>
#include <chargefw/charges/charge_set.h>
#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/io/charge_result.h>
#include <chargefw/io/mol2_input.h>
#include <chargefw/io/mol2_output.h>
#include <chargefw/version.h>
#include <snitch/snitch.hpp>

#include <limits>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

namespace io = chargefw::io;
namespace charges = chargefw::charges;
namespace core = chargefw::core;

using chargefw::test::calculation_result;

auto generated_record() -> io::ImportedMoleculeRecord {
    return {.molecule =
                core::Molecule{{core::Atom{6, 0, "C1"}, core::Atom{8, -1, "O1"}},
                               {core::Bond{0, 1, core::BondOrder::DOUBLE}},
                               {core::Conformer{{core::Position{0.0, 0.0, 0.0},
                                                 core::Position{1.23456789012345, 0.0, 0.0}}},
                                core::Conformer{{core::Position{0.1, 0.0, 0.0},
                                                 core::Position{1.33456789012345, 0.0, 0.0}}}},
                               "carbonyl"},
            .identity = {.source_name = "input", .record_id = "record"}};
}

[[nodiscard]] auto occurrences(const std::string_view text, const std::string_view value)
    -> std::size_t {
    std::size_t result = 0;
    for (auto position = text.find(value); position != std::string_view::npos;
         position = text.find(value, position + value.size())) {
        ++result;
    }
    return result;
}

// The atom types written for a molecule, in atom order.
[[nodiscard]] auto written_atom_types(std::vector<core::Atom> atoms, std::vector<core::Bond> bonds)
    -> std::vector<std::string> {
    const auto atom_count = atoms.size();
    auto record = io::ImportedMoleculeRecord{
        .molecule = core::Molecule{std::move(atoms),
                                   std::move(bonds),
                                   {core::Conformer{std::vector<core::Position>(atom_count)}}},
        .identity = {.source_name = "input"}};
    const auto result = calculation_result(
        {std::move(record)},
        charges::ChargeSet{"formal",
                           {{.target = {.molecule_index = 0},
                             .charges = charges::AtomicCharges{std::vector(atom_count, 0.0)}}}});
    auto output = std::ostringstream{};
    io::Mol2Writer{output}.write(result);

    auto input = std::istringstream{output.str()};
    auto types = std::vector<std::string>{};
    auto line = std::string{};
    while (std::getline(input, line) && line != "@<TRIPOS>ATOM") {
    }
    for (std::size_t index = 0; index < atom_count && std::getline(input, line); ++index) {
        auto fields = std::istringstream{line};
        auto field = std::string{};
        for (auto column = 0; column < 6; ++column) {
            fields >> field;
        }
        types.push_back(field);
    }
    return types;
}

} // namespace

TEST_CASE("MOL2 output writes SYBYL atom types derived from the topology", "[io][mol2]") {
    using core::Atom;
    using core::Bond;
    constexpr auto single = core::BondOrder::SINGLE;
    constexpr auto double_ = core::BondOrder::DOUBLE;
    constexpr auto triple = core::BondOrder::TRIPLE;
    struct Case {
        std::string_view name;
        std::vector<Atom> atoms;
        std::vector<Bond> bonds;
        std::vector<std::string> types;
    };
    const auto cases = std::vector<Case>{
        {.name = "acetate shares the charge between its oxygens",
         .atoms = {Atom{6}, Atom{6}, Atom{8}, Atom{8, -1}, Atom{1}, Atom{1}, Atom{1}},
         .bonds = {Bond{0, 1, single}, Bond{1, 2, double_}, Bond{1, 3, single}, Bond{0, 4, single},
                   Bond{0, 5, single}, Bond{0, 6, single}},
         .types = {"C.3", "C.2", "O.co2", "O.co2", "H", "H", "H"}},
        {.name = "acetic acid keeps its carbonyl and hydroxyl oxygens",
         .atoms = {Atom{6}, Atom{6}, Atom{8}, Atom{8}, Atom{1}, Atom{1}, Atom{1}, Atom{1}},
         .bonds = {Bond{0, 1, single}, Bond{1, 2, double_}, Bond{1, 3, single}, Bond{0, 4, single},
                   Bond{0, 5, single}, Bond{0, 6, single}, Bond{3, 7, single}},
         .types = {"C.3", "C.2", "O.2", "O.3", "H", "H", "H", "H"}},
        {.name = "carbonate is not a carboxylate",
         .atoms = {Atom{6}, Atom{8}, Atom{8, -1}, Atom{8, -1}},
         .bonds = {Bond{0, 1, double_}, Bond{0, 2, single}, Bond{0, 3, single}},
         .types = {"C.2", "O.2", "O.3", "O.3"}},
        {.name = "methyl phosphate dianion",
         .atoms = {Atom{15}, Atom{8}, Atom{8, -1}, Atom{8, -1}, Atom{8}, Atom{6}, Atom{1}, Atom{1},
                   Atom{1}},
         .bonds = {Bond{0, 1, double_}, Bond{0, 2, single}, Bond{0, 3, single}, Bond{0, 4, single},
                   Bond{4, 5, single}, Bond{5, 6, single}, Bond{5, 7, single}, Bond{5, 8, single}},
         .types = {"P.3", "O.co2", "O.co2", "O.co2", "O.3", "C.3", "H", "H", "H"}},
        {.name = "methylammonium",
         .atoms = {Atom{6}, Atom{7, 1}, Atom{1}, Atom{1}, Atom{1}, Atom{1}, Atom{1}, Atom{1}},
         .bonds = {Bond{0, 1, single}, Bond{0, 2, single}, Bond{0, 3, single}, Bond{0, 4, single},
                   Bond{1, 5, single}, Bond{1, 6, single}, Bond{1, 7, single}},
         .types = {"C.3", "N.4", "H", "H", "H", "H", "H", "H"}},
        {.name = "methylamine",
         .atoms = {Atom{6}, Atom{7}, Atom{1}, Atom{1}, Atom{1}, Atom{1}, Atom{1}},
         .bonds = {Bond{0, 1, single}, Bond{0, 2, single}, Bond{0, 3, single}, Bond{0, 4, single},
                   Bond{1, 5, single}, Bond{1, 6, single}},
         .types = {"C.3", "N.3", "H", "H", "H", "H", "H"}},
        {.name = "N-methylacetamide",
         .atoms = {Atom{6}, Atom{6}, Atom{8}, Atom{7}, Atom{6}, Atom{1}, Atom{1}, Atom{1}, Atom{1},
                   Atom{1}, Atom{1}, Atom{1}},
         .bonds = {Bond{0, 1, single}, Bond{1, 2, double_}, Bond{1, 3, single}, Bond{3, 4, single},
                   Bond{0, 5, single}, Bond{0, 6, single}, Bond{0, 7, single}, Bond{3, 8, single},
                   Bond{4, 9, single}, Bond{4, 10, single}, Bond{4, 11, single}},
         .types = {"C.3", "C.2", "O.2", "N.am", "C.3", "H", "H", "H", "H", "H", "H", "H"}},
        {.name = "vinylamine is conjugated",
         .atoms = {Atom{6}, Atom{6}, Atom{7}, Atom{1}, Atom{1}, Atom{1}, Atom{1}, Atom{1}},
         .bonds = {Bond{0, 1, double_}, Bond{1, 2, single}, Bond{0, 3, single}, Bond{0, 4, single},
                   Bond{1, 5, single}, Bond{2, 6, single}, Bond{2, 7, single}},
         .types = {"C.2", "C.2", "N.pl3", "H", "H", "H", "H", "H"}},
        {.name = "benzene is written in its Kekule form",
         .atoms = {Atom{6}, Atom{6}, Atom{6}, Atom{6}, Atom{6}, Atom{6}, Atom{1}, Atom{1}, Atom{1},
                   Atom{1}, Atom{1}, Atom{1}},
         .bonds = {Bond{0, 1, double_}, Bond{1, 2, single}, Bond{2, 3, double_}, Bond{3, 4, single},
                   Bond{4, 5, double_}, Bond{5, 0, single}, Bond{0, 6, single}, Bond{1, 7, single},
                   Bond{2, 8, single}, Bond{3, 9, single}, Bond{4, 10, single},
                   Bond{5, 11, single}},
         .types = {"C.2", "C.2", "C.2", "C.2", "C.2", "C.2", "H", "H", "H", "H", "H", "H"}},
        {.name = "nitromethane",
         .atoms = {Atom{6}, Atom{7, 1}, Atom{8}, Atom{8, -1}, Atom{1}, Atom{1}, Atom{1}},
         .bonds = {Bond{0, 1, single}, Bond{1, 2, double_}, Bond{1, 3, single}, Bond{0, 4, single},
                   Bond{0, 5, single}, Bond{0, 6, single}},
         .types = {"C.3", "N.pl3", "O.2", "O.3", "H", "H", "H"}},
        {.name = "acetonitrile",
         .atoms = {Atom{6}, Atom{6}, Atom{7}, Atom{1}, Atom{1}, Atom{1}},
         .bonds = {Bond{0, 1, single}, Bond{1, 2, triple}, Bond{0, 3, single}, Bond{0, 4, single},
                   Bond{0, 5, single}},
         .types = {"C.3", "C.1", "N.1", "H", "H", "H"}},
        {.name = "guanidinium",
         .atoms = {Atom{6}, Atom{7, 1}, Atom{7}, Atom{7}, Atom{1}, Atom{1}, Atom{1}, Atom{1},
                   Atom{1}, Atom{1}},
         .bonds = {Bond{0, 1, double_}, Bond{0, 2, single}, Bond{0, 3, single}, Bond{1, 4, single},
                   Bond{1, 5, single}, Bond{2, 6, single}, Bond{2, 7, single}, Bond{3, 8, single},
                   Bond{3, 9, single}},
         .types = {"C.cat", "N.pl3", "N.pl3", "N.pl3", "H", "H", "H", "H", "H", "H"}},
        {.name = "dimethyl sulfoxide",
         .atoms = {Atom{16}, Atom{8}, Atom{6}, Atom{6}, Atom{1}, Atom{1}, Atom{1}, Atom{1}, Atom{1},
                   Atom{1}},
         .bonds = {Bond{0, 1, double_}, Bond{0, 2, single}, Bond{0, 3, single}, Bond{2, 4, single},
                   Bond{2, 5, single}, Bond{2, 6, single}, Bond{3, 7, single}, Bond{3, 8, single},
                   Bond{3, 9, single}},
         .types = {"S.O", "O.2", "C.3", "C.3", "H", "H", "H", "H", "H", "H"}},
        {.name = "dimethyl sulfoxide with separated charges",
         .atoms = {Atom{16, 1}, Atom{8, -1}, Atom{6}, Atom{6}, Atom{1}, Atom{1}, Atom{1}, Atom{1},
                   Atom{1}, Atom{1}},
         .bonds = {Bond{0, 1, single}, Bond{0, 2, single}, Bond{0, 3, single}, Bond{2, 4, single},
                   Bond{2, 5, single}, Bond{2, 6, single}, Bond{3, 7, single}, Bond{3, 8, single},
                   Bond{3, 9, single}},
         .types = {"S.O", "O.3", "C.3", "C.3", "H", "H", "H", "H", "H", "H"}},
        {.name = "dimethyl sulfone",
         .atoms = {Atom{16}, Atom{8}, Atom{8}, Atom{6}, Atom{6}, Atom{1}, Atom{1}, Atom{1}, Atom{1},
                   Atom{1}, Atom{1}},
         .bonds = {Bond{0, 1, double_}, Bond{0, 2, double_}, Bond{0, 3, single}, Bond{0, 4, single},
                   Bond{3, 5, single}, Bond{3, 6, single}, Bond{3, 7, single}, Bond{4, 8, single},
                   Bond{4, 9, single}, Bond{4, 10, single}},
         .types = {"S.O2", "O.2", "O.2", "C.3", "C.3", "H", "H", "H", "H", "H", "H"}},
        {.name = "methanesulfenate is not a sulfoxide",
         .atoms = {Atom{6}, Atom{16}, Atom{8, -1}, Atom{1}, Atom{1}, Atom{1}},
         .bonds = {Bond{0, 1, single}, Bond{1, 2, single}, Bond{0, 3, single}, Bond{0, 4, single},
                   Bond{0, 5, single}},
         .types = {"C.3", "S.3", "O.3", "H", "H", "H"}},
        {.name = "sodium chloride",
         .atoms = {Atom{11, 1}, Atom{17, -1}},
         .bonds = {},
         .types = {"Na", "Cl"}},
    };

    for (const auto& item : cases) {
        CAPTURE(item.name);
        CHECK(written_atom_types(item.atoms, item.bonds) == item.types);
    }
}

TEST_CASE("result-owned MOL2 writes every conformer with its charges", "[io][mol2]") {
    const auto result = calculation_result(
        {generated_record()},
        charges::ChargeSet{
            "qeq",
            {{.target = {.molecule_index = 0, .conformer_index = 0},
              .charges = charges::AtomicCharges{{0.123456789012345, -0.123456789012345}}},
             {.target = {.molecule_index = 0, .conformer_index = 1},
              .charges = charges::AtomicCharges{{0.25, -0.25}}}},
            "qeq-default"});
    auto output = std::ostringstream{};
    io::Mol2Writer{output}.write(result);
    const auto text = output.str();

    CHECK(occurrences(text, "@<TRIPOS>MOLECULE") == 2);
    CHECK(text.contains("record_conformer_1"));
    CHECK(text.contains("record_conformer_2"));
    CHECK(text.contains("1.23456789012345"));
    CHECK(text.contains(" 0.1 0 0 "));
    CHECK(text.contains("Generated by ChargeFW " CHARGEFW_VERSION_STRING));
    CHECK(text.contains("; method=qeq; parameter_set=qeq-default"));
    CHECK(text.contains("@<TRIPOS>BOND\n1 1 2 2"));
    CHECK(occurrences(text, "@<TRIPOS>SUBSTRUCTURE\n1 UNL 1") == 2);
    const auto second_record = text.find("@<TRIPOS>MOLECULE", 1);
    REQUIRE(second_record != std::string::npos);
    const auto first_text = std::string_view{text}.substr(0, second_record);
    const auto second_text = std::string_view{text}.substr(second_record);
    CHECK(first_text.contains(" 0.1235\n"));
    CHECK_FALSE(first_text.contains(" 0.2500\n"));
    CHECK(second_text.contains(" 0.2500\n"));
    CHECK_FALSE(second_text.contains(" 0.1235\n"));

    auto round_trip_stream = std::istringstream{text};
    auto reader = io::Mol2Reader{round_trip_stream};
    const auto first = reader.next();
    const auto second = reader.next();
    REQUIRE(first.has_value());
    REQUIRE(second.has_value());
    CHECK(first->identity.record_id == "record_conformer_1");
    CHECK(second->identity.record_id == "record_conformer_2");
    CHECK(first->molecule.atom_count() == 2);
    CHECK(first->molecule.atom(1).formal_charge() == -1);
    CHECK(second->molecule.bond_count() == 1);
    CHECK(first->molecule.conformer(0)[1].x == 1.23456789012345);
    CHECK(second->molecule.conformer(0)[1].x == 1.33456789012345);
    CHECK_FALSE(reader.next().has_value());
}

TEST_CASE("result-owned MOL2 broadcasts molecule charges to retained conformers", "[io][mol2]") {
    auto second = generated_record();
    second.identity.record_id = "second";
    const auto result = calculation_result(
        {generated_record(), std::move(second)},
        charges::ChargeSet{
            "formal",
            {{.target = {.molecule_index = 0}, .charges = charges::AtomicCharges{{0.25, -0.25}}},
             {.target = {.molecule_index = 1}, .charges = charges::AtomicCharges{{0.5, -0.5}}}}});
    auto output = std::ostringstream{};
    io::Mol2Writer{output}.write(result);

    CHECK(occurrences(output.str(), " 0.2500\n") == 2);
    CHECK(occurrences(output.str(), " -0.2500\n") == 2);
    CHECK(occurrences(output.str(), " 0.5000\n") == 2);
    CHECK(occurrences(output.str(), " -0.5000\n") == 2);
}

TEST_CASE("result-owned MOL2 rejects failed results and invalid values before output",
          "[io][mol2]") {
    {
        const auto cancelled = io::make_charge_calculation_result(
            {generated_record()}, {},
            {.status = chargefw::calculation::ExecutionStatus::cancelled});
        auto output = std::ostringstream{};
        CHECK_THROWS_AS(io::Mol2Writer{output}.write(cancelled), std::invalid_argument);
        CHECK(output.str().empty());
    }

    {
        auto record = generated_record();
        record.molecule = core::Molecule{{core::Atom{6}, core::Atom{8}}};
        const auto missing = calculation_result(
            {std::move(record)},
            charges::ChargeSet{"formal",
                               {{.target = {.molecule_index = 0},
                                 .charges = charges::AtomicCharges{{0.0, 0.0}}}}});
        auto output = std::ostringstream{};
        CHECK_THROWS_AS(io::Mol2Writer{output}.write(missing), std::invalid_argument);
        CHECK(output.str().empty());
    }

    {
        auto record = generated_record();
        record.molecule = core::Molecule{{}, {}, {core::Conformer{{}}}};
        const auto empty = calculation_result(
            {std::move(record)}, charges::ChargeSet{"formal",
                                                    {{.target = {.molecule_index = 0},
                                                      .charges = charges::AtomicCharges{{}}}}});
        auto output = std::ostringstream{};
        CHECK_THROWS_AS(io::Mol2Writer{output}.write(empty), std::invalid_argument);
        CHECK(output.str().empty());
    }

    {
        auto record = generated_record();
        record.molecule = core::Molecule{
            {core::Atom{6}, core::Atom{8}},
            {},
            {core::Conformer{{{.x = std::numeric_limits<double>::infinity(), .y = 0.0, .z = 0.0},
                              {.x = 1.0, .y = 0.0, .z = 0.0}}}}};
        const auto nonfinite_coordinates = calculation_result(
            {std::move(record)},
            charges::ChargeSet{"formal",
                               {{.target = {.molecule_index = 0, .conformer_index = 0},
                                 .charges = charges::AtomicCharges{{0.0, 0.0}}}}});
        auto output = std::ostringstream{};
        CHECK_THROWS_AS(io::Mol2Writer{output}.write(nonfinite_coordinates), std::invalid_argument);
        CHECK(output.str().empty());
    }
}
