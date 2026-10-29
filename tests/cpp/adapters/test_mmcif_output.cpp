#include "support/test_charge_results.h"

#include <chargefw/adapters/charge_result.h>
#include <chargefw/adapters/gemmi/mmcif_input.h>
#include <chargefw/adapters/gemmi/mmcif_output.h>
#include <chargefw/adapters/gemmi/pdb_input.h>
#include <chargefw/charges/atomic_charges.h>
#include <chargefw/charges/charge_collection.h>

#include <gemmi/cif.hpp>
#include <gemmi/read_cif.hpp>
#include <snitch/snitch.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

namespace {

namespace adapters = chargefw::adapters;
namespace charges = chargefw::charges;
namespace core = chargefw::core;
namespace mmcif_input = chargefw::adapters::gemmi::mmcif_input;
namespace mmcif_output = chargefw::adapters::gemmi::mmcif_output;
namespace pdb_input = chargefw::adapters::gemmi::pdb_input;

using chargefw::test::calculation_result;

auto generated_record(std::string id = "molecule") -> adapters::ImportedMoleculeRecord {
    return {
        .molecule =
            core::Molecule{
                {core::Atom{6}, core::Atom{8, -1}},
                {core::Bond{0, 1, core::BondOrder::DOUBLE}},
                {core::Conformer{{core::Position{0.0, 0.0, 0.0}, core::Position{1.2, 0.0, 0.0}}},
                 core::Conformer{{core::Position{0.1, 0.0, 0.0}, core::Position{1.3, 0.0, 0.0}}}}},
        .identity = {.source = "input", .record_id = std::move(id)}};
}

auto structural_result() -> adapters::ChargeCalculationResult {
    auto record = generated_record("structural");
    auto first_sites = std::vector<adapters::SourceAtomReference>{
        {.position = 2,
         .id = "Csite",
         .structural_labels =
             adapters::SourceStructuralLabels{
                 .author = {.atom = "AUTH_C", .residue = "ARES", .chain = "AC", .sequence = "17"},
                 .label = {.atom = "LABEL_C", .residue = "LRES", .chain = "LC", .sequence = "7"},
                 .entity = "E1",
                 .alternate_location = "A"}},
        {.position = 3,
         .id = "001",
         .structural_labels = adapters::SourceStructuralLabels{
             .author = {.atom = "AUTH_O", .residue = "ARES", .chain = "AC", .sequence = "17"},
             .label = {.atom = "LABEL_O", .residue = "LRES", .chain = "LC", .sequence = "7"},
             .entity = "E1"}}};
    auto second_sites = first_sites;
    second_sites[0].id = "9223372036854775808";
    second_sites[0].structural_labels->alternate_location = "B";
    second_sites[1].id = "other";
    record.import_metadata = adapters::MoleculeImportMetadata{
        .format = adapters::MolecularSourceFormat::mmcif,
        .atoms = first_sites,
        .conformers = {{.position = 0, .id = "source-model-1", .sites = std::move(first_sites)},
                       {.position = 1, .id = "source-model-2", .sites = std::move(second_sites)}},
        .record_selection = "all",
        .alternate_location_selection = "first-source-order",
        .conformer_selection = "all",
        .bond_strategy = "none"};
    return calculation_result(
        {std::move(record)},
        charges::ChargeSet{
            "qeq",
            {{.target = {.molecule_index = 0, .conformer_index = 0},
              .charges = charges::AtomicCharges{{0.123456789012345, -0.123456789012345}}},
             {.target = {.molecule_index = 0, .conformer_index = 1},
              .charges = charges::AtomicCharges{{0.2, -0.2}}}},
            "qeq-default"});
}

struct LabelledAtom {
    std::string residue;
    std::string sequence;
    std::string atom;
    std::string author_chain = "A";
    std::string author_atom; // the label atom name when empty
};

// Carbon atoms with mmCIF labels in label chain A, split into the given residue instances.
auto labelled_record(const std::vector<LabelledAtom>& atoms, std::vector<core::Bond> bonds,
                     std::vector<adapters::SourceComponentInstance> components)
    -> adapters::ImportedMoleculeRecord {
    auto positions = std::vector<core::Position>{};
    auto sites = std::vector<adapters::SourceAtomReference>{};
    for (std::size_t index = 0; index < atoms.size(); ++index) {
        positions.push_back({.x = static_cast<double>(index), .y = 0.0, .z = 0.0});
        const auto& atom = atoms[index];
        const auto label = adapters::SourceHierarchyLabels{
            .atom = atom.atom, .residue = atom.residue, .chain = "A", .sequence = atom.sequence};
        auto author = label;
        author.chain = atom.author_chain;
        if (!atom.author_atom.empty()) {
            author.atom = atom.author_atom;
        }
        sites.push_back({.position = index,
                         .structural_labels =
                             adapters::SourceStructuralLabels{.author = author, .label = label}});
    }
    auto record = adapters::ImportedMoleculeRecord{
        .molecule = core::Molecule{std::vector<core::Atom>(atoms.size(), core::Atom{6}),
                                   std::move(bonds),
                                   {core::Conformer{std::move(positions)}}},
        .identity = {.source = "input", .record_id = "labelled"}};
    record.import_metadata =
        adapters::MoleculeImportMetadata{.format = adapters::MolecularSourceFormat::mmcif,
                                         .atoms = sites,
                                         .conformers = {{.position = 0, .sites = std::move(sites)}},
                                         .components = std::move(components)};
    return record;
}

[[nodiscard]] auto formal(const std::size_t atom_count) -> charges::ChargeSet {
    return charges::ChargeSet{
        "formal",
        {{.target = {.molecule_index = 0},
          .charges = charges::AtomicCharges{std::vector<double>(atom_count, 0.0)}}}};
}

[[nodiscard]] auto bond_set(const core::Molecule& molecule)
    -> std::set<std::tuple<std::size_t, std::size_t, core::BondOrder>> {
    auto result = std::set<std::tuple<std::size_t, std::size_t, core::BondOrder>>{};
    for (const auto& bond : molecule.bonds()) {
        result.emplace(std::min(bond.first_atom_index(), bond.second_atom_index()),
                       std::max(bond.first_atom_index(), bond.second_atom_index()), bond.order());
    }
    return result;
}

} // namespace

TEST_CASE("result-owned mmCIF preserves hierarchy and joins conformer charges",
          "[adapters][mmcif]") {
    const auto result = structural_result();
    auto output = std::ostringstream{};
    mmcif_output::MmcifWriter{output}.write(result, "test");
    auto document = ::gemmi::cif::read_string(output.str());
    REQUIRE(document.blocks.size() == 1);
    auto& block = document.blocks.front();
    auto component_bonds =
        block.find("_chem_comp_bond.", {"comp_id", "atom_id_1", "atom_id_2", "value_order"});
    REQUIRE(component_bonds.length() == 1);
    CHECK(::gemmi::cif::as_string(component_bonds[0][0]) == "LRES");
    CHECK(::gemmi::cif::as_string(component_bonds[0][1]) == "LABEL_C");
    CHECK(::gemmi::cif::as_string(component_bonds[0][2]) == "LABEL_O");
    CHECK(::gemmi::cif::as_string(component_bonds[0][3]) == "doub");
    CHECK_FALSE(block.has_mmcif_category("_struct_conn."));
    // _chem_comp would need component types that ChargeFW does not know, so it is omitted.
    CHECK_FALSE(block.has_mmcif_category("_chem_comp."));

    auto sites = block.find("_atom_site.", {"id", "label_atom_id", "label_comp_id", "label_asym_id",
                                            "auth_atom_id", "auth_comp_id", "auth_asym_id",
                                            "label_alt_id", "pdbx_PDB_model_num"});
    REQUIRE(sites.length() == 4);
    CHECK(::gemmi::cif::as_string(sites[0][0]) == "1");
    CHECK(::gemmi::cif::as_string(sites[3][0]) == "4");
    CHECK(::gemmi::cif::as_string(sites[0][1]) == "LABEL_C");
    CHECK(::gemmi::cif::as_string(sites[0][3]) == "LC");
    CHECK(::gemmi::cif::as_string(sites[0][4]) == "AUTH_C");
    CHECK(::gemmi::cif::as_string(sites[0][6]) == "AC");
    CHECK(::gemmi::cif::as_string(sites[0][7]) == "A");
    CHECK(::gemmi::cif::as_string(sites[2][7]) == "B");
    CHECK(::gemmi::cif::as_string(sites[2][8]) == "2");
    CHECK_FALSE(output.str().contains("Csite"));
    CHECK_FALSE(output.str().contains("9223372036854775808"));

    auto values = block.find("_sb_ncbr_partial_atomic_charges.", {"type_id", "atom_id", "charge"});
    REQUIRE(values.length() == 4);
    auto metadata =
        block.find("_sb_ncbr_partial_atomic_charges_meta.", {"software_name", "software_version"});
    // One charge set covers both conformers under one id.
    REQUIRE(metadata.length() == 1);
    CHECK(::gemmi::cif::as_string(metadata[0][0]) == "ChargeFW");
    CHECK(::gemmi::cif::as_string(metadata[0][1]) == "test");
    CHECK(::gemmi::cif::as_string(values[0][0]) == "1");
    CHECK(::gemmi::cif::as_string(values[2][0]) == "1");
    CHECK(::gemmi::cif::as_string(values[2][1]) == "3");
    CHECK(::gemmi::cif::as_string(values[0][2]) == "0.1235");

    auto round_trip_stream = std::istringstream{output.str()};
    auto round_trip_reader = mmcif_input::MmcifReader{
        round_trip_stream, {}, {.bond_strategy = adapters::gemmi::BondStrategy::explicit_bonds}};
    const auto round_trip = round_trip_reader.next();
    REQUIRE(round_trip.has_value());
    CHECK(round_trip->molecule.atom_count() == 2);
    CHECK(round_trip->molecule.conformer_count() == 2);
    CHECK(bond_set(round_trip->molecule) == bond_set(result.inputs().front().molecule));
}

TEST_CASE("mmCIF bonds between PDB residues read back unchanged", "[adapters][mmcif]") {
    // A peptide link from the templates and a disulfide from CONECT. PDB residues have no label
    // sequence numbers, so the links must be found by author chain and number.
    auto pdb = std::istringstream{
        "ATOM      1  N   CYS A   1      -0.966   0.493   1.500  1.00  0.00           N\n"
        "ATOM      2  CA  CYS A   1       0.257   0.418   0.692  1.00  0.00           C\n"
        "ATOM      3  C   CYS A   1      -0.094   0.017  -0.716  1.00  0.00           C\n"
        "ATOM      4  O   CYS A   1      -1.056  -0.682  -0.923  1.00  0.00           O\n"
        "ATOM      5  CB  CYS A   1       1.204  -0.620   1.296  1.00  0.00           C\n"
        "ATOM      6  SG  CYS A   1       2.100  -1.300   0.000  1.00  0.00           S\n"
        "ATOM      7  N   CYS A   2       0.661   0.439  -1.742  1.00  0.00           N\n"
        "ATOM      8  CA  CYS A   2       0.343   0.091  -3.135  1.00  0.00           C\n"
        "ATOM      9  C   CYS A   2       1.418   0.656  -4.054  1.00  0.00           C\n"
        "ATOM     10  O   CYS A   2       2.398   1.233  -3.587  1.00  0.00           O\n"
        "ATOM     11  CB  CYS A   2      -1.034   0.627  -3.509  1.00  0.00           C\n"
        "ATOM     12  SG  CYS A   2       2.000  -1.200  -2.000  1.00  0.00           S\n"
        "CONECT    6   12\n"
        "CONECT   12    6\n"
        "END\n"};
    auto record = pdb_input::PdbReader{pdb,
                                       "input.pdb",
                                       {.bond_strategy = adapters::gemmi::BondStrategy::hybrid}}
                      .next();
    REQUIRE(record.has_value());
    const auto expected = bond_set(record->molecule);
    const auto atom_count = record->molecule.atom_count();
    auto output = std::ostringstream{};
    mmcif_output::MmcifWriter{output}.write(
        calculation_result({std::move(*record)}, formal(atom_count)));

    auto document = ::gemmi::cif::read_string(output.str());
    auto connections = document.blocks.front().find("_struct_conn.", {"conn_type_id"});
    auto types = std::set<std::string>{};
    for (const auto row : connections) {
        types.insert(::gemmi::cif::as_string(row[0]));
    }
    CHECK(types == std::set<std::string>{"covale", "disulf"});

    auto input = std::istringstream{output.str()};
    const auto round_trip =
        mmcif_input::MmcifReader{
            input, {}, {.bond_strategy = adapters::gemmi::BondStrategy::explicit_bonds}}
            .next();
    REQUIRE(round_trip.has_value());
    CHECK(bond_set(round_trip->molecule) == expected);
}

TEST_CASE("mmCIF bond output rejects bonds that residue definitions cannot express",
          "[adapters][mmcif]") {
    const auto two_ligands = std::vector<LabelledAtom>{
        {"LIG", "1", "C"}, {"LIG", "1", "O"}, {"LIG", "2", "C"}, {"LIG", "2", "O"}};
    const auto separate =
        std::vector<adapters::SourceComponentInstance>{{"LIG", {0, 1}}, {"LIG", {2, 3}}};
    const auto double_bonds = std::vector<core::Bond>{core::Bond{0, 1, core::BondOrder::DOUBLE},
                                                      core::Bond{2, 3, core::BondOrder::DOUBLE}};
    // Residues are told apart by the reader's residue list, here by author chain alone.
    const auto other_author_chain = std::vector<LabelledAtom>{
        {"LIG", "1", "C"}, {"LIG", "1", "O"}, {"LIG", "1", "C", "B"}, {"LIG", "1", "O", "B"}};
    for (const auto& atoms : {two_ligands, other_author_chain}) {
        auto output = std::ostringstream{};
        mmcif_output::MmcifWriter{output}.write(
            calculation_result({labelled_record(atoms, double_bonds, separate)}, formal(4)));
        auto document = ::gemmi::cif::read_string(output.str());
        CHECK(document.blocks.front().find("_chem_comp_bond.", {"comp_id"}).length() == 1);
    }

    struct Case {
        std::string_view description;
        adapters::ImportedMoleculeRecord record;
    };
    const auto cases = std::vector<Case>{
        {"repeated atom name", labelled_record({{"LIG", "1", "C"}, {"LIG", "1", "C"}},
                                               {core::Bond{0, 1}}, {{"LIG", {0, 1}}})},
        {"different orders", labelled_record(two_ligands,
                                             {core::Bond{0, 1, core::BondOrder::DOUBLE},
                                              core::Bond{2, 3, core::BondOrder::SINGLE}},
                                             separate)},
        {"bond in one residue only", labelled_record(two_ligands, {core::Bond{0, 1}}, separate)},
        {"bond between residues to an atom with different author and label names",
         labelled_record({{"ALA", "1", "C", "A", "AUTH_C"}, {"ALA", "2", "N"}}, {core::Bond{0, 1}},
                         {{"ALA", {0}}, {"ALA", {1}}})},
    };
    for (const auto& [description, record] : cases) {
        CAPTURE(description);
        const auto result = calculation_result({record}, formal(record.molecule.atom_count()));
        auto output = std::ostringstream{};
        CHECK_THROWS_AS(mmcif_output::MmcifWriter{output}.write(result), std::invalid_argument);
        CHECK(output.str().empty());
    }
}

TEST_CASE("mmCIF holds one charge set per result for the same molecules", "[adapters][mmcif]") {
    const auto charge_set = [](std::string method, std::optional<std::string> parameter_set,
                               const double charge) {
        return charges::ChargeSet{std::move(method),
                                  {{.target = {.molecule_index = 0},
                                    .charges = charges::AtomicCharges{{charge, -charge}}}},
                                  std::move(parameter_set)};
    };
    const auto qeq = calculation_result({generated_record()}, charge_set("qeq", "QEq", 0.25));
    const auto eem = calculation_result({generated_record()}, charge_set("eem", "EEM", 0.5));
    auto output = std::ostringstream{};
    mmcif_output::MmcifWriter{output}.write(std::array{&qeq, &eem}, "test");
    auto document = ::gemmi::cif::read_string(output.str());
    auto& block = document.blocks.front();
    auto metadata =
        block.find("_sb_ncbr_partial_atomic_charges_meta.", {"id", "method", "parameter_set"});
    REQUIRE(metadata.length() == 2);
    CHECK(::gemmi::cif::as_string(metadata[1][0]) == "2");
    CHECK(::gemmi::cif::as_string(metadata[1][1]) == "eem");
    CHECK(::gemmi::cif::as_string(metadata[1][2]) == "EEM");
    // Two atoms in two conformers for each set, each set with its own charges.
    auto charges_by_site = std::map<std::pair<std::string, std::string>, std::string>{};
    for (const auto row :
         block.find("_sb_ncbr_partial_atomic_charges.", {"type_id", "atom_id", "charge"})) {
        charges_by_site[{::gemmi::cif::as_string(row[0]), ::gemmi::cif::as_string(row[1])}] =
            ::gemmi::cif::as_string(row[2]);
    }
    CHECK(charges_by_site.size() == 8);
    CHECK(charges_by_site[{"1", "3"}] == "0.2500");
    CHECK(charges_by_site[{"2", "3"}] == "0.5000");
    CHECK(charges_by_site[{"2", "4"}] == "-0.5000");

    const auto other_molecule =
        calculation_result({generated_record("other")}, charge_set("eem", "EEM", 0.5));
    const auto same_set = calculation_result({generated_record()}, charge_set("qeq", "QEq", 0.5));
    using Results = std::array<const adapters::ChargeCalculationResult*, 2>;
    for (const auto& rejected :
         {Results{&qeq, &other_molecule}, Results{&qeq, &same_set}, Results{&qeq, nullptr}}) {
        auto rejected_output = std::ostringstream{};
        CHECK_THROWS_AS(mmcif_output::MmcifWriter{rejected_output}.write(rejected),
                        std::invalid_argument);
        CHECK(rejected_output.str().empty());
    }
}

TEST_CASE("result-owned mmCIF broadcasts molecule charges and uses fresh unique blocks",
          "[adapters][mmcif]") {
    auto first = generated_record("same");
    auto second = first;
    const auto result = calculation_result(
        {std::move(first), std::move(second)},
        charges::ChargeSet{
            "formal",
            {{.target = {.molecule_index = 0}, .charges = charges::AtomicCharges{{0.25, -0.25}}},
             {.target = {.molecule_index = 1}, .charges = charges::AtomicCharges{{0.1, -0.1}}}}});
    auto output = std::ostringstream{};
    mmcif_output::MmcifWriter{output}.write(result);
    auto document = ::gemmi::cif::read_string(output.str());
    REQUIRE(document.blocks.size() == 2);
    CHECK(document.blocks[0].name == "same");
    CHECK(document.blocks[1].name == "same_2");
    auto metadata = document.blocks[0].find("_sb_ncbr_partial_atomic_charges_meta.", {"id"});
    REQUIRE(metadata.length() == 1);
    auto charge_rows =
        document.blocks[0].find("_sb_ncbr_partial_atomic_charges.", {"atom_id", "charge"});
    REQUIRE(charge_rows.length() == 4);
    CHECK(::gemmi::cif::as_string(charge_rows[0][0]) == "1");
    CHECK(::gemmi::cif::as_string(charge_rows[2][0]) == "3");
    CHECK(::gemmi::cif::as_string(charge_rows[0][1]) == ::gemmi::cif::as_string(charge_rows[2][1]));
}

TEST_CASE("result-owned mmCIF rejects failed results geometry and out-of-range charges",
          "[adapters][mmcif]") {
    const auto cancelled = adapters::make_charge_calculation_result(
        {generated_record("cancelled")}, {},
        {.status = chargefw::calculation::ExecutionStatus::cancelled});
    auto output = std::ostringstream{};
    CHECK_THROWS_AS(mmcif_output::MmcifWriter{output}.write(cancelled), std::invalid_argument);

    auto coordinate_free = generated_record("missing");
    coordinate_free.molecule = core::Molecule{{core::Atom{6}, core::Atom{8}}};
    const auto missing_result =
        calculation_result({std::move(coordinate_free)},
                           charges::ChargeSet{"formal",
                                              {{.target = {.molecule_index = 0},
                                                .charges = charges::AtomicCharges{{0.0, 0.0}}}}});
    CHECK_THROWS_AS(mmcif_output::MmcifWriter{output}.write(missing_result), std::invalid_argument);

    const auto range_result =
        calculation_result({generated_record("range")},
                           charges::ChargeSet{"formal",
                                              {{.target = {.molecule_index = 0},
                                                .charges = charges::AtomicCharges{{5.1, -5.1}}}}});
    CHECK_THROWS_AS(mmcif_output::MmcifWriter{output}.write(range_result), std::invalid_argument);

    auto nonfinite = generated_record("nonfinite");
    nonfinite.molecule = core::Molecule{
        {core::Atom{6}, core::Atom{8}},
        {},
        {core::Conformer{{{.x = std::numeric_limits<double>::infinity(), .y = 0.0, .z = 0.0},
                          {.x = 1.0, .y = 0.0, .z = 0.0}}}}};
    const auto nonfinite_result = calculation_result(
        {std::move(nonfinite)},
        charges::ChargeSet{"formal",
                           {{.target = {.molecule_index = 0, .conformer_index = 0},
                             .charges = charges::AtomicCharges{{0.0, 0.0}}}}});
    CHECK_THROWS_AS(mmcif_output::MmcifWriter{output}.write(nonfinite_result),
                    std::invalid_argument);
    CHECK(output.str().empty());
}
