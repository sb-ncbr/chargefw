#include "support/test_charge_results.h"

#include <chargefw/charges/atomic_charges.h>
#include <chargefw/charges/charge_set.h>
#include <chargefw/io/charge_result.h>
#include <chargefw/io/mmcif_input.h>
#include <chargefw/io/mmcif_output.h>
#include <chargefw/io/pdb_input.h>
#include <chargefw/version.h>

#include <gemmi/cif.hpp>
#include <gemmi/read_cif.hpp>
#include <snitch/snitch.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <filesystem>
#include <format>
#include <fstream>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <ranges>
#include <set>
#include <sstream>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

namespace {

namespace io = chargefw::io;
namespace charges = chargefw::charges;
namespace core = chargefw::core;

using chargefw::test::calculation_result;

auto generated_record(std::string id = "molecule") -> io::ImportedMoleculeRecord {
    return {
        .molecule =
            core::Molecule{
                {core::Atom{6}, core::Atom{8, -1}},
                {core::Bond{0, 1, core::BondOrder::DOUBLE}},
                {core::Conformer{{core::Position{0.0, 0.0, 0.0}, core::Position{1.2, 0.0, 0.0}}},
                 core::Conformer{{core::Position{0.1, 0.0, 0.0}, core::Position{1.3, 0.0, 0.0}}}}},
        .identity = {.source_name = "input", .record_id = std::move(id)}};
}

auto structural_result() -> io::ChargeCalculationResult {
    auto record = generated_record("structural");
    auto first_sites = std::vector<io::SourceAtomReference>{
        {.position = 2,
         .id = "Csite",
         .structural_labels =
             std::make_shared<const io::SourceStructuralLabels>(io::SourceStructuralLabels{
                 .author = {.atom = "AUTH_C", .residue = "ARES", .chain = "AC", .sequence = "17"},
                 .label = {.atom = "LABEL_C", .residue = "LRES", .chain = "LC", .sequence = "7"},
                 .entity = "E1",
                 .alternate_location = "A"})},
        {.position = 3,
         .id = "001",
         .structural_labels =
             std::make_shared<const io::SourceStructuralLabels>(io::SourceStructuralLabels{
                 .author = {.atom = "AUTH_O", .residue = "ARES", .chain = "AC", .sequence = "17"},
                 .label = {.atom = "LABEL_O", .residue = "LRES", .chain = "LC", .sequence = "7"},
                 .entity = "E1"})}};
    auto second_sites = first_sites;
    second_sites[0].id = "9223372036854775808";
    auto second_labels = *second_sites[0].structural_labels;
    second_labels.alternate_location = "B";
    second_sites[0].structural_labels =
        std::make_shared<const io::SourceStructuralLabels>(std::move(second_labels));
    second_sites[1].id = "other";
    record.import_metadata = io::MoleculeImportMetadata{
        .format = io::MolecularSourceFormat::mmcif,
        .atoms = first_sites,
        .conformers = {{.position = 0, .id = "source-model-1", .sites = std::move(first_sites)},
                       {.position = 1, .id = "source-model-2", .sites = std::move(second_sites)}},
        .policy = {.selection = "all",
                   .alternate_locations = "first-source-order",
                   .conformers = "all",
                   .bonds = "none"}};
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
                     std::vector<io::SourceComponentInstance> components)
    -> io::ImportedMoleculeRecord {
    auto positions = std::vector<core::Position>{};
    auto sites = std::vector<io::SourceAtomReference>{};
    for (std::size_t index = 0; index < atoms.size(); ++index) {
        positions.push_back({.x = static_cast<double>(index), .y = 0.0, .z = 0.0});
        const auto& atom = atoms[index];
        const auto label = io::SourceHierarchyLabels{
            .atom = atom.atom, .residue = atom.residue, .chain = "A", .sequence = atom.sequence};
        auto author = label;
        author.chain = atom.author_chain;
        if (!atom.author_atom.empty()) {
            author.atom = atom.author_atom;
        }
        sites.push_back({.position = index,
                         .structural_labels = std::make_shared<const io::SourceStructuralLabels>(
                             io::SourceStructuralLabels{.author = author, .label = label})});
    }
    auto record = io::ImportedMoleculeRecord{
        .molecule = core::Molecule{std::vector<core::Atom>(atoms.size(), core::Atom{6}),
                                   std::move(bonds),
                                   {core::Conformer{std::move(positions)}}},
        .identity = {.source_name = "input", .record_id = "labelled"}};
    record.import_metadata =
        io::MoleculeImportMetadata{.format = io::MolecularSourceFormat::mmcif,
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

TEST_CASE("result-owned mmCIF preserves hierarchy and joins conformer charges", "[io][mmcif]") {
    const auto result = structural_result();
    auto output = std::ostringstream{};
    io::MmcifWriter{output}.write(result);
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
    CHECK(::gemmi::cif::as_string(metadata[0][1]) == CHARGEFW_VERSION_STRING);
    CHECK(::gemmi::cif::as_string(values[0][0]) == "1");
    CHECK(::gemmi::cif::as_string(values[2][0]) == "1");
    CHECK(::gemmi::cif::as_string(values[2][1]) == "3");
    CHECK(::gemmi::cif::as_string(values[0][2]) == "0.1235");

    auto round_trip_stream = std::istringstream{output.str()};
    auto round_trip_reader =
        io::MmcifReader{round_trip_stream, {}, {.bonds = io::BondStrategy::explicit_bonds}};
    const auto round_trip = round_trip_reader.next();
    REQUIRE(round_trip.has_value());
    CHECK(round_trip->molecule.atom_count() == 2);
    CHECK(round_trip->molecule.conformer_count() == 2);
    CHECK(bond_set(round_trip->molecule) == bond_set(result.inputs().front().molecule));
}

TEST_CASE("mmCIF bonds between PDB residues read back unchanged", "[io][mmcif]") {
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
    auto record = io::PdbReader{pdb, "input.pdb", {.bonds = io::BondStrategy::hybrid}}.next();
    REQUIRE(record.has_value());
    const auto expected = bond_set(record->molecule);
    const auto atom_count = record->molecule.atom_count();
    auto output = std::ostringstream{};
    io::MmcifWriter{output}.write(calculation_result({std::move(*record)}, formal(atom_count)));

    auto document = ::gemmi::cif::read_string(output.str());
    auto connections = document.blocks.front().find("_struct_conn.", {"conn_type_id"});
    auto types = std::set<std::string>{};
    for (const auto row : connections) {
        types.insert(::gemmi::cif::as_string(row[0]));
    }
    CHECK(types == std::set<std::string>{"covale", "disulf"});

    auto input = std::istringstream{output.str()};
    const auto round_trip =
        io::MmcifReader{input, {}, {.bonds = io::BondStrategy::explicit_bonds}}.next();
    REQUIRE(round_trip.has_value());
    CHECK(bond_set(round_trip->molecule) == expected);
}

TEST_CASE("mmCIF bonds of mixed protonation states read back and write again unchanged",
          "[io][mmcif]") {
    struct Residue {
        std::string_view name;
        std::string_view atoms;
        std::map<std::string_view, int> charges;
    };
    // Histidines protonated on ND1, on NE2, and on both with +1 on NE2, and arginines with and
    // without HH12, so that residues of one name have their double bonds on different atoms.
    const auto residues = std::vector<Residue>{
        {"HIS", "CB CG ND1 CD2 CE1 NE2 HD2 HE1 HD1", {}},
        {"HIS", "CB CG ND1 CD2 CE1 NE2 HD2 HE1 HE2", {}},
        {"HIS", "CB CG ND1 CD2 CE1 NE2 HD2 HE1 HD1 HE2", {{"NE2", 1}}},
        {"ARG", "CG CD HD2 HD3 NE CZ NH1 NH2 HE HH11 HH12 HH21 HH22", {{"NH2", 1}}},
        {"ARG", "CG CD HD2 HD3 NE CZ NH1 NH2 HE HH11 HH21 HH22", {}},
    };
    auto pdb = std::string{};
    auto serial = 0;
    for (std::size_t residue_index = 0; residue_index < residues.size(); ++residue_index) {
        const auto& [name, atoms, charges] = residues[residue_index];
        for (const auto atom_range : std::views::split(atoms, ' ')) {
            const auto atom = std::string{std::string_view{atom_range}};
            const auto charge = charges.contains(atom) ? charges.at(atom) : 0;
            ++serial;
            pdb += std::format("ATOM  {:>5} {:<4} {} A{:>4}    {:8.3f}{:8.3f}   0.000  1.00 20.00"
                               "          {:>2}{}\n",
                               serial, atom.size() == 4 ? atom : " " + atom, name,
                               2 * residue_index + 1, static_cast<double>(serial),
                               static_cast<double>(residue_index), atom.substr(0, 1),
                               charge == 0 ? "  " : "1+");
        }
    }
    pdb += "END\n";

    const auto write = [](io::ImportedMoleculeRecord record) {
        const auto atom_count = record.molecule.atom_count();
        auto output = std::ostringstream{};
        io::MmcifWriter{output}.write(calculation_result({std::move(record)}, formal(atom_count)));
        return output.str();
    };
    const auto read_mmcif = [](const std::string& text) {
        auto input = std::istringstream{text};
        auto record = io::MmcifReader{input, {}, {.bonds = io::BondStrategy::hybrid}}.next();
        REQUIRE(record.has_value());
        return std::move(*record);
    };

    auto pdb_input = std::istringstream{pdb};
    const auto original = io::PdbReader{pdb_input, {}, {.bonds = io::BondStrategy::hybrid}}.next();
    REQUIRE(original.has_value());
    const auto written = write(*original);

    auto delocalized = std::set<std::string>{};
    auto document = ::gemmi::cif::read_string(written);
    for (const auto row : document.blocks.front().find(
             "_chem_comp_bond.", {"comp_id", "atom_id_1", "atom_id_2", "value_order"})) {
        if (row.str(3) == "delo") {
            const auto [first, second] = std::minmax({row.str(1), row.str(2)});
            delocalized.insert(std::format("{} {}-{}", row.str(0), first, second));
        }
    }
    CHECK(delocalized ==
          std::set<std::string>{"HIS CE1-ND1", "HIS CE1-NE2", "ARG CZ-NH1", "ARG CZ-NH2"});

    const auto round_trip = read_mmcif(written);
    CHECK(bond_set(round_trip.molecule) == bond_set(original->molecule));
    for (std::size_t index = 0; index < original->molecule.atom_count(); ++index) {
        CHECK(round_trip.molecule.atom(index).formal_charge() ==
              original->molecule.atom(index).formal_charge());
    }
    CHECK(write(round_trip) == written);
}

TEST_CASE("mmCIF bond output rejects bonds that residue definitions cannot express",
          "[io][mmcif]") {
    const auto two_ligands = std::vector<LabelledAtom>{
        {"LIG", "1", "C"}, {"LIG", "1", "O"}, {"LIG", "2", "C"}, {"LIG", "2", "O"}};
    const auto separate =
        std::vector<io::SourceComponentInstance>{{"LIG", {0, 1}}, {"LIG", {2, 3}}};
    const auto double_bonds = std::vector<core::Bond>{core::Bond{0, 1, core::BondOrder::DOUBLE},
                                                      core::Bond{2, 3, core::BondOrder::DOUBLE}};
    // Residues are told apart by the reader's residue list, here by author chain alone.
    const auto other_author_chain = std::vector<LabelledAtom>{
        {"LIG", "1", "C"}, {"LIG", "1", "O"}, {"LIG", "1", "C", "B"}, {"LIG", "1", "O", "B"}};
    for (const auto& atoms : {two_ligands, other_author_chain}) {
        auto output = std::ostringstream{};
        io::MmcifWriter{output}.write(
            calculation_result({labelled_record(atoms, double_bonds, separate)}, formal(4)));
        auto document = ::gemmi::cif::read_string(output.str());
        CHECK(document.blocks.front().find("_chem_comp_bond.", {"comp_id"}).length() == 1);
    }

    struct Case {
        std::string_view description;
        io::ImportedMoleculeRecord record;
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
        CHECK_THROWS_AS(io::MmcifWriter{output}.write(result), std::invalid_argument);
        CHECK(output.str().empty());
    }
}

TEST_CASE("mmCIF holds one charge set per result for the same molecules", "[io][mmcif]") {
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
    io::MmcifWriter{output}.write(std::array{&qeq, &eem});
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
    using Results = std::array<const io::ChargeCalculationResult*, 2>;
    for (const auto& rejected :
         {Results{&qeq, &other_molecule}, Results{&qeq, &same_set}, Results{&qeq, nullptr}}) {
        auto rejected_output = std::ostringstream{};
        CHECK_THROWS_AS(io::MmcifWriter{rejected_output}.write(rejected), std::invalid_argument);
        CHECK(rejected_output.str().empty());
    }
}

TEST_CASE("mmCIF combines results for one structure read twice", "[io][mmcif]") {
    // Each reading creates its own structural label objects; records compare by label values.
    const auto read = [] {
        auto pdb = std::istringstream{
            "ATOM      1  O   HOH A   1       0.000   0.000   0.000  1.00  0.00           O\n"
            "ATOM      2  H1  HOH A   1       0.957   0.000   0.000  1.00  0.00           H\n"
            "END\n"};
        auto record = io::PdbReader{pdb, "input.pdb"}.next();
        REQUIRE(record.has_value());
        return std::move(*record);
    };
    const auto charge_set = [](std::string method, const double charge) {
        return charges::ChargeSet{std::move(method),
                                  {{.target = {.molecule_index = 0},
                                    .charges = charges::AtomicCharges{{-charge, charge}}}}};
    };
    const auto qeq = calculation_result({read()}, charge_set("qeq", 0.4));
    const auto eem = calculation_result({read()}, charge_set("eem", 0.3));
    auto output = std::ostringstream{};
    CHECK_NOTHROW(io::MmcifWriter{output}.write(std::array{&qeq, &eem}));

    auto renamed = read();
    auto labels = *renamed.import_metadata->atoms[1].structural_labels;
    labels.author.atom = "H2";
    renamed.import_metadata->atoms[1].structural_labels =
        std::make_shared<const io::SourceStructuralLabels>(std::move(labels));
    const auto other = calculation_result({std::move(renamed)}, charge_set("eem", 0.3));
    auto rejected_output = std::ostringstream{};
    CHECK_THROWS_AS(io::MmcifWriter{rejected_output}.write(std::array{&qeq, &other}),
                    std::invalid_argument);
}

TEST_CASE("mmCIF annotation adds charges for the selected sites to the source structure",
          "[io][mmcif]") {
    // Gapped site IDs, an alternate location, water, a second model, and a block without atoms.
    // Only the charged sites need integer IDs; the water of the second model has a text ID.
    constexpr auto source = std::string_view{R"(data_notes
_note.text 'kept'
data_structure
_cell.length_a 10.0
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.pdbx_PDB_model_num
HETATM 10 C C1 . LIG A 1 . 0.0 0.0 0.0 1
HETATM 20 O O1 A LIG A 1 . 1.2 0.0 0.0 1
HETATM 21 O O1 B LIG A 1 . 0.0 1.2 0.0 1
HETATM 30 O O  . HOH B 2 . 5.0 5.0 5.0 1
HETATM 40 C C1 . LIG A 1 . 0.1 0.0 0.0 2
HETATM 50 O O1 A LIG A 1 . 1.3 0.0 0.0 2
HETATM 51 O O1 B LIG A 1 . 0.1 1.2 0.0 2
HETATM waterO O O  . HOH B 2 . 5.1 5.0 5.0 2
)"};
    const auto read = [](const std::string_view text) {
        auto input = std::istringstream{std::string{text}};
        // The reader is gone before the charges are written.
        auto record = io::MmcifReader{input,
                                      "input.cif",
                                      {.selection = io::RecordSelection::polymers_and_ligands,
                                       .conformers = io::ConformerSelection::first}}
                          .next();
        REQUIRE(record.has_value());
        return std::move(*record);
    };
    const auto charge_set = [](std::string method) {
        return charges::ChargeSet{
            std::move(method),
            {{.target = {.molecule_index = 0}, .charges = charges::AtomicCharges{{0.25, -0.25}}}}};
    };
    const auto result = calculation_result({read(source)}, charge_set("qeq"));
    auto output = std::ostringstream{};
    io::MmcifWriter{output}.annotate(std::array{&result});

    auto document = gemmi::cif::read_string(output.str());
    REQUIRE(document.blocks.size() == 2);
    CHECK(*document.blocks[0].find_value("_note.text") == "'kept'");
    auto& block = document.blocks[1];
    CHECK(*block.find_value("_cell.length_a") == "10.0");
    auto site_ids = std::vector<std::string>{};
    for (const auto row : block.find("_atom_site.", {"id"})) {
        site_ids.push_back(row.str(0));
    }
    CHECK(site_ids == std::vector<std::string>{"10", "20", "21", "30", "40", "50", "51", "waterO"});
    auto charged = std::vector<std::pair<std::string, std::string>>{};
    for (const auto row :
         block.find("_sb_ncbr_partial_atomic_charges.", {"type_id", "atom_id", "charge"})) {
        CHECK(row.str(0) == "1");
        charged.emplace_back(row.str(1), row.str(2));
    }
    CHECK(charged ==
          std::vector<std::pair<std::string, std::string>>{{"10", "0.2500"}, {"20", "-0.2500"}});

    // The output is input again: the same molecule, so the same settings give the same charges.
    const auto reread = read(output.str());
    CHECK(reread.molecule == result.inputs().front().molecule);
    CHECK(reread.import_metadata == result.inputs().front().import_metadata);

    const auto again = calculation_result({reread}, charge_set("eem"));
    auto rejected = std::ostringstream{};
    CHECK_THROWS_AS(io::MmcifWriter{rejected}.annotate(std::array{&again}), std::invalid_argument);
    auto replaced = std::ostringstream{};
    io::MmcifWriter{replaced}.annotate(std::array{&again}, io::ExistingCharges::replace);
    auto replaced_document = gemmi::cif::read_string(replaced.str());
    auto methods =
        replaced_document.blocks[1].find("_sb_ncbr_partial_atomic_charges_meta.", {"method"});
    REQUIRE(methods.length() == 1);
    CHECK(methods[0].str(0) == "eem");

    // Appending keeps the set that is there and continues its IDs; the same method is refused.
    auto appended = std::ostringstream{};
    io::MmcifWriter{appended}.annotate(std::array{&again}, io::ExistingCharges::append);
    auto appended_document = gemmi::cif::read_string(appended.str());
    auto sets = std::vector<std::pair<std::string, std::string>>{};
    for (const auto row : appended_document.blocks[1].find("_sb_ncbr_partial_atomic_charges_meta.",
                                                           {"id", "method"})) {
        sets.emplace_back(row.str(0), row.str(1));
    }
    CHECK(sets == std::vector<std::pair<std::string, std::string>>{{"1", "qeq"}, {"2", "eem"}});
    auto set_of_row = std::vector<std::string>{};
    for (const auto row :
         appended_document.blocks[1].find("_sb_ncbr_partial_atomic_charges.", {"type_id"})) {
        set_of_row.push_back(row.str(0));
    }
    CHECK(set_of_row == std::vector<std::string>{"1", "1", "2", "2"});
    const auto repeated = calculation_result({reread}, charge_set("qeq"));
    CHECK_THROWS_AS(
        io::MmcifWriter{rejected}.annotate(std::array{&repeated}, io::ExistingCharges::append),
        std::invalid_argument);
    // Without charges in the source, appending only adds; charges of dictionary 1.0 are refused.
    auto first = std::ostringstream{};
    io::MmcifWriter{first}.annotate(std::array{&result}, io::ExistingCharges::append);
    CHECK(first.str() == output.str());
    const auto legacy = calculation_result(
        {read(std::string{source} + "_sb_ncbr_partial_atomic_charges_meta.id 1\n"
                                    "_sb_ncbr_partial_atomic_charges_meta.type empirical\n"
                                    "_sb_ncbr_partial_atomic_charges_meta.method 'qeq/set'\n")},
        charge_set("eem"));
    CHECK_THROWS_AS(
        io::MmcifWriter{rejected}.annotate(std::array{&legacy}, io::ExistingCharges::append),
        std::invalid_argument);

    // Existing sets may lack the optional parameter set and have any text as ID; the new set
    // follows the largest number among the IDs.
    const auto with_charges = [&](const std::string_view declaration,
                                  const std::string_view set_rows, const std::string_view rows) {
        return std::format("{}{}loop_\n_sb_ncbr_partial_atomic_charges_meta.id\n"
                           "_sb_ncbr_partial_atomic_charges_meta.type\n"
                           "_sb_ncbr_partial_atomic_charges_meta.method\n"
                           "_sb_ncbr_partial_atomic_charges_meta.software_name\n"
                           "_sb_ncbr_partial_atomic_charges_meta.software_version\n{}"
                           "loop_\n_sb_ncbr_partial_atomic_charges.type_id\n"
                           "_sb_ncbr_partial_atomic_charges.atom_id\n"
                           "_sb_ncbr_partial_atomic_charges.charge\n{}",
                           source, declaration, set_rows, rows);
    };
    const auto append_to = [&](const std::string& text) {
        const auto added = calculation_result({read(text)}, charge_set("eem"));
        auto written = std::ostringstream{};
        io::MmcifWriter{written}.annotate(std::array{&added}, io::ExistingCharges::append);
        return written.str();
    };
    constexpr auto two_sets = std::string_view{"qm_ref empirical qm other 1\n"
                                               "2147483647 empirical qeq other 1\n"};
    auto text_ids = gemmi::cif::read_string(
        append_to(with_charges("", two_sets, "qm_ref 10 0.1\n2147483647 30 0.2\n")));
    auto id_rows = std::vector<std::pair<std::string, std::string>>{};
    for (const auto row : text_ids.blocks[1].find("_sb_ncbr_partial_atomic_charges_meta.",
                                                  {"id", "method", "parameter_set"})) {
        id_rows.emplace_back(row.str(0), row.str(1));
        CHECK_FALSE(row.has2(2));
    }
    CHECK(id_rows == std::vector<std::pair<std::string, std::string>>{
                         {"qm_ref", "qm"}, {"2147483647", "qeq"}, {"2147483648", "eem"}});

    // Refused: a charge of an unknown site, two charges of one site in a set, and tables declared
    // to follow another version of the dictionary.
    constexpr auto v10 = std::string_view{"_audit_conform.dict_name mmcif_charges_v10.dic\n"};
    for (const auto& text : {with_charges("", two_sets, "qm_ref 99999 0.1\n"),
                             with_charges("", two_sets, "qm_ref 10 0.1\nqm_ref 10 0.2\n"),
                             with_charges(v10, two_sets, "qm_ref 10 0.1\n")}) {
        CHECK_THROWS_AS(append_to(text), std::invalid_argument);
    }

    // Equal molecules from inputs that differ elsewhere are not one source, in either order.
    auto other_source = std::string{source};
    other_source.replace(other_source.find("10.0"), 4, "20.0");
    const auto other = calculation_result({read(other_source)}, charge_set("eem"));
    REQUIRE(other.inputs().front() == result.inputs().front());
    CHECK_THROWS_AS(io::MmcifWriter{rejected}.annotate(std::array{&result, &other}),
                    std::invalid_argument);
    CHECK_THROWS_AS(io::MmcifWriter{rejected}.annotate(std::array{&other, &result}),
                    std::invalid_argument);
    // Separate reads of the same text are.
    const auto same = calculation_result({read(source)}, charge_set("eem"));
    auto combined = std::ostringstream{};
    CHECK_NOTHROW(io::MmcifWriter{combined}.annotate(std::array{&result, &same}));

    // A declaration of another dictionary as single items stays, without a location.
    auto declared_source = std::string{source};
    declared_source.insert(declared_source.find("_cell."),
                           "_audit_conform.dict_name mmcif_pdbx.dic\n"
                           "_audit_conform.dict_version 5.0\n");
    const auto declared = calculation_result({read(declared_source)}, charge_set("eem"));
    auto declared_output = std::ostringstream{};
    io::MmcifWriter{declared_output}.annotate(std::array{&declared});
    auto declared_document = gemmi::cif::read_string(declared_output.str());
    auto declarations = std::vector<std::pair<std::string, std::string>>{};
    for (const auto row :
         declared_document.blocks[1].find("_audit_conform.", {"dict_name", "dict_version"})) {
        declarations.emplace_back(row.str(0), row.str(1));
    }
    CHECK(declarations == std::vector<std::pair<std::string, std::string>>{
                              {"mmcif_pdbx.dic", "5.0"}, {"mmcif_charges_v11.dic", "1.1"}});

    // A molecule that was not read from mmCIF has no structure to annotate.
    const auto generated = calculation_result({generated_record()}, formal(2));
    CHECK_THROWS_AS(io::MmcifWriter{rejected}.annotate(std::array{&generated}),
                    std::invalid_argument);
}

TEST_CASE("PDB annotation puts the charges on the calculated atoms of the converted structure",
          "[io][mmcif][pdb]") {
    // The block name, the number of atom sites and, in calculation order, the sites of the charge
    // rows, after checking that each is the site of its calculated atom.
    struct Annotated {
        std::string block_name;
        std::size_t site_count;
        std::vector<std::string> charged_ids;
    };
    const auto annotate = [](const std::string& text, const io::RecordSelection selection) {
        auto pdb = std::istringstream{text};
        auto read = io::PdbReader{pdb, "input.pdb", {.selection = selection}}.next();
        REQUIRE(read.has_value());
        const auto atom_count = read->molecule.atom_count();
        const auto result = calculation_result({std::move(*read)}, formal(atom_count));
        const auto& molecule = result.inputs().front().molecule;
        auto output = std::ostringstream{};
        io::MmcifWriter{output}.annotate(std::array{&result});

        auto document = gemmi::cif::read_string(output.str());
        auto& block = document.sole_block();
        auto positions = std::map<std::string, core::Position>{};
        for (const auto row : block.find("_atom_site.", {"id", "Cartn_x", "Cartn_y", "Cartn_z"})) {
            positions.emplace(row.str(0), core::Position{std::stod(row[1]), std::stod(row[2]),
                                                         std::stod(row[3])});
        }
        auto charges = block.find("_sb_ncbr_partial_atomic_charges.", {"atom_id"});
        REQUIRE(charges.length() == atom_count);
        auto charged_ids = std::vector<std::string>{};
        for (std::size_t atom = 0; atom < atom_count; ++atom) {
            charged_ids.push_back(charges[static_cast<int>(atom)].str(0));
            CHECK(positions.at(charged_ids.back()) == molecule.conformer(0)[atom]);
        }
        return Annotated{.block_name = block.name,
                         .site_count = positions.size(),
                         .charged_ids = std::move(charged_ids)};
    };

    // Two alanines, two ligands, and water: the whole structure is written, and the selected atoms
    // get charges.
    const auto fixture = [] {
        auto file = std::ifstream{std::filesystem::path{CHARGEFW_TEST_SOURCE_DIR} / "tests" /
                                  "fixtures" / "synthetic" / "pdb" / "annotation.pdb"};
        REQUIRE(file.good());
        return std::string{std::istreambuf_iterator<char>{file}, std::istreambuf_iterator<char>{}};
    }();
    const auto polymers = annotate(fixture, io::RecordSelection::polymers);
    CHECK(polymers.site_count == 16);
    CHECK(polymers.charged_ids.size() == 11);
    // The block is named as the molecule, by the HEADER ID code, not by the file.
    CHECK(polymers.block_name == "1ANN");

    // Gemmi joins the two parts of the first ligand, so its sites are not in the order of the
    // lines. The atom that differs from another only in the segment is a site with a charge too.
    const auto reordered = std::string{
        R"(HETATM    1  C1  LIG A   1       0.000   0.000   0.000  1.00 20.00      SEG1 C
HETATM    2  O   HOH A   2       5.000   0.000   0.000  1.00 20.00      SEG1 O
HETATM    3  O1  LIG A   1       1.200   0.000   0.000  1.00 20.00      SEG1 O
HETATM    4  C1  LIG A   1       0.000   5.000   0.000  1.00 20.00      SEG2 C
END
)"};
    CHECK(annotate(reordered, io::RecordSelection::all).charged_ids ==
          std::vector<std::string>{"1", "3", "2", "4"});
}

TEST_CASE("result-owned mmCIF broadcasts molecule charges and uses fresh unique blocks",
          "[io][mmcif]") {
    auto first = generated_record("same");
    auto second = first;
    const auto result = calculation_result(
        {std::move(first), std::move(second)},
        charges::ChargeSet{
            "formal",
            {{.target = {.molecule_index = 0}, .charges = charges::AtomicCharges{{0.25, -0.25}}},
             {.target = {.molecule_index = 1}, .charges = charges::AtomicCharges{{0.1, -0.1}}}}});
    auto output = std::ostringstream{};
    io::MmcifWriter{output}.write(result);
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

    // A suffix is skipped when a record already has that name.
    auto records = std::vector<io::ImportedMoleculeRecord>{};
    auto assignments = std::vector<charges::ChargeAssignment>{};
    for (const auto* name : {"same", "same_2", "same", "same"}) {
        assignments.push_back({.target = {.molecule_index = records.size()},
                               .charges = charges::AtomicCharges{{0.0, 0.0}}});
        records.push_back(generated_record(name));
    }
    auto named_output = std::ostringstream{};
    io::MmcifWriter{named_output}.write(calculation_result(
        std::move(records), charges::ChargeSet{"formal", std::move(assignments)}));
    auto named = ::gemmi::cif::read_string(named_output.str());
    auto names = std::vector<std::string>{};
    for (const auto& block : named.blocks) {
        names.push_back(block.name);
    }
    CHECK(names == std::vector<std::string>{"same", "same_2", "same_3", "same_4"});
}

TEST_CASE("result-owned mmCIF rejects failed results geometry and out-of-range charges",
          "[io][mmcif]") {
    const auto cancelled = io::make_charge_calculation_result(
        {generated_record("cancelled")}, {},
        {.status = chargefw::calculation::ExecutionStatus::cancelled});
    auto output = std::ostringstream{};
    CHECK_THROWS_AS(io::MmcifWriter{output}.write(cancelled), std::invalid_argument);

    auto coordinate_free = generated_record("missing");
    coordinate_free.molecule = core::Molecule{{core::Atom{6}, core::Atom{8}}};
    const auto missing_result =
        calculation_result({std::move(coordinate_free)},
                           charges::ChargeSet{"formal",
                                              {{.target = {.molecule_index = 0},
                                                .charges = charges::AtomicCharges{{0.0, 0.0}}}}});
    CHECK_THROWS_AS(io::MmcifWriter{output}.write(missing_result), std::invalid_argument);

    // The dictionary range is inclusive.
    const auto range_result = [](const double charge) {
        return calculation_result(
            {generated_record("range")},
            charges::ChargeSet{"formal",
                               {{.target = {.molecule_index = 0},
                                 .charges = charges::AtomicCharges{{charge, -charge}}}}});
    };
    auto in_range = std::ostringstream{};
    CHECK_NOTHROW(io::MmcifWriter{in_range}.write(range_result(10.0)));
    CHECK_THROWS_AS(io::MmcifWriter{output}.write(range_result(10.1)), std::invalid_argument);

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
    CHECK_THROWS_AS(io::MmcifWriter{output}.write(nonfinite_result), std::invalid_argument);
    CHECK(output.str().empty());
}
