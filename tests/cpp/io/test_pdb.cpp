#include "support/test_exceptions.h"

#include <chargefw/core/bond.h>
#include <chargefw/io/component_dictionary.h>
#include <chargefw/io/conformer_selection.h>
#include <chargefw/io/molecule_record.h>
#include <chargefw/io/parse_error.h>
#include <chargefw/io/pdb_input.h>
#include <chargefw/io/structure_options.h>
#include <optional>
#include <snitch/snitch.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <exception>
#include <format>
#include <map>
#include <memory>
#include <ranges>
#include <set>
#include <sstream>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace io = chargefw::io;

namespace {

[[nodiscard]] auto read_first(const std::string& text, const io::StructureOptions& options = {}) {
    std::istringstream input{text};
    auto reader = io::PdbReader{input, {}, options};
    auto record = reader.next();
    REQUIRE(record.has_value());
    return std::move(*record);
}

// PDB atom lines of one residue from space-separated atom names, whose first letter is the element.
[[nodiscard]] auto pdb_residue(const std::string_view residue, const std::string_view atoms,
                               const std::map<std::string_view, int>& charges = {}) {
    auto text = std::string{};
    auto serial = 0;
    for (const auto name_range : std::views::split(atoms, ' ')) {
        const auto name = std::string{std::string_view{name_range}};
        const auto charge = charges.contains(name) ? charges.at(name) : 0;
        const auto charge_field =
            charge == 0 ? "  " : std::format("{}{}", std::abs(charge), charge > 0 ? '+' : '-');
        ++serial;
        text += std::format("ATOM  {:>5} {:<4} {:>3} A   1    {:8.3f}   0.000   0.000  1.00 20.00"
                            "          {:>2}{}\n",
                            serial, name.size() == 4 ? name : " " + name, residue,
                            static_cast<double>(serial), name.substr(0, 1), charge_field);
    }
    return text + "END\n";
}

} // namespace

static_assert(!std::is_copy_constructible_v<io::PdbReader> &&
              !std::is_copy_assignable_v<io::PdbReader>);
static_assert(std::is_move_constructible_v<io::PdbReader>);

TEST_CASE("PDB templates supplement N-terminal hydrogen connectivity", "[io][pdb]") {
    using Strategy = io::BondStrategy;
    for (const auto strategy : {Strategy::templates, Strategy::hybrid}) {
        std::istringstream input{
            R"pdb(ATOM      1  N   ALA A   7       0.000   0.000   0.000  1.00 20.00           N1+
ATOM      2  CA  ALA A   7       1.000   0.000   0.000  1.00 20.00           C
ATOM      3  H1  ALA A   7       0.000   1.000   0.000  1.00 20.00           H
ATOM      4  H2  ALA A   7       0.000   0.000   1.000  1.00 20.00           H
ATOM      5  H3  ALA A   7       0.000  -1.000   0.000  1.00 20.00           H
CONECT    1    3
END
)pdb"};
        auto reader = io::PdbReader{input, {}, {.bonds = strategy}};
        const auto record = reader.next();
        REQUIRE(record.has_value());
        CHECK(record->molecule.atom_count() == 5);
        CHECK(record->molecule.bond_count() == 4);
        CHECK(record->molecule.atom(0).formal_charge() == 1);
        for (const auto& bond : record->molecule.bonds()) {
            CHECK(bond.first_atom_index() == 0);
            CHECK(bond.order() == chargefw::core::BondOrder::SINGLE);
        }
    }

    // Hydride names H1 the hydrogen it puts on the carboxyl carbon of a C-terminal residue without
    // OXT. Only the first residue of the peptide is an N-terminus, so H1 is not bonded to its N.
    const auto record = read_first(
        R"pdb(ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00 20.00           N
ATOM      2  CA  GLY A   1       1.450   0.000   0.000  1.00 20.00           C
ATOM      3  C   GLY A   1       2.000   1.400   0.000  1.00 20.00           C
ATOM      4  O   GLY A   1       1.300   2.400   0.000  1.00 20.00           O
ATOM      5  N   GLY A   2       3.300   1.500   0.000  1.00 20.00           N
ATOM      6  CA  GLY A   2       3.900   2.800   0.000  1.00 20.00           C
ATOM      7  C   GLY A   2       5.400   2.800   0.000  1.00 20.00           C
ATOM      8  O   GLY A   2       6.000   1.700   0.000  1.00 20.00           O
ATOM      9  H1  GLY A   2       5.900   3.800   0.000  1.00 20.00           H
END
)pdb",
        {.bonds = Strategy::templates});
    for (const auto& bond : record.molecule.bonds()) {
        CHECK(bond.first_atom_index() != 8);
        CHECK(bond.second_atom_index() != 8);
    }

    // Without the C of the first residue there is no peptide bond, so the second residue is an
    // N-terminus and its H1 is bonded to its N.
    const auto unlinked = read_first(
        R"pdb(ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00 20.00           N
ATOM      2  CA  GLY A   1       1.450   0.000   0.000  1.00 20.00           C
ATOM      3  N   GLY A   2       3.300   1.500   0.000  1.00 20.00           N
ATOM      4  CA  GLY A   2       3.900   2.800   0.000  1.00 20.00           C
ATOM      5  H1  GLY A   2       2.900   0.600   0.000  1.00 20.00           H
END
)pdb",
        {.bonds = Strategy::templates});
    CHECK(std::ranges::any_of(unlinked.molecule.bonds(), [](const auto& bond) {
        return bond.first_atom_index() == 2 && bond.second_atom_index() == 4;
    }));
}

TEST_CASE("PDB templates follow the protonation state of histidine and arginine", "[io][pdb]") {
    struct Case {
        std::string_view residue;
        std::string_view atoms;
        std::map<std::string_view, int> charges;
        std::set<std::string> double_bonds;
    };
    // The charged forms follow the wwPDB protonation variants, such as HIS_LL_DHE2 for histidine
    // without HE2. Without formal charges in the input, only neutral forms fit.
    const auto cases = std::vector<Case>{
        {"HIS", "CB CG ND1 CD2 CE1 NE2 HD2 HE1 HD1", {}, {"CG-CD2", "CE1-NE2"}},
        {"HIS",
         "CB CG ND1 CD2 CE1 NE2 HD2 HE1 HD1",
         {{"ND1", 1}, {"NE2", -1}},
         {"CG-CD2", "ND1-CE1"}},
        {"HIS", "CB CG ND1 CD2 CE1 NE2 HD2 HE1 HE2", {}, {"CG-CD2", "ND1-CE1"}},
        {"HIS", "CB CG ND1 CD2 CE1 NE2 HD2 HE1 HD1 HE2", {{"ND1", 1}}, {"CG-CD2", "ND1-CE1"}},
        {"HIS", "CB CG ND1 CD2 CE1 NE2 HD2 HE1 HD1 HE2", {{"NE2", 1}}, {"CG-CD2", "CE1-NE2"}},
        {"ARG", "CG CD HD2 HD3 NE CZ NH1 NH2 HE HH11 HH12 HH21 HH22", {{"NH2", 1}}, {"CZ-NH2"}},
        {"ARG", "CG CD HD2 HD3 NE CZ NH1 NH2 HE HH11 HH12 HH21 HH22", {{"NH1", 1}}, {"CZ-NH1"}},
        {"ARG", "CG CD HD2 HD3 NE CZ NH1 NH2 HE HH11 HH21 HH22", {}, {"CZ-NH1"}},
        {"ARG", "CG CD HD2 HD3 NE CZ NH1 NH2 HE HH11 HH12 HH21", {}, {"CZ-NH2"}},
        {"ARG", "CG CD HD2 HD3 NE CZ NH1 NH2 HH11 HH12 HH21 HH22", {}, {"NE-CZ"}},
        // No structure fits a cation without its formal charge or a group missing a bonded atom.
        {"HIS", "CB CG ND1 CD2 CE1 NE2 HD2 HE1 HD1 HE2", {}, {"CG-CD2", "ND1-CE1"}},
        {"HIS", "CG ND1 CD2 CE1 NE2 HD2 HE1 HD1", {}, {"CG-CD2", "ND1-CE1"}},
        {"ARG", "CG CD HD2 HD3 NE CZ NH1 NH2 HE HH11 HH12 HH21 HH22", {}, {"CZ-NH2"}},
        {"ARG", "NE CZ NH1 NH2 HE HH11 HH12 HH21 HH22", {}, {"CZ-NH2"}},
        {"ARG", "CG CD HD2 NE CZ NH1 NH2 HH11 HH12 HH21 HH22", {}, {"CZ-NH2"}},
    };

    for (const auto& [residue, atoms, charges, double_bonds] : cases) {
        const auto text = pdb_residue(residue, atoms, charges);
        const auto record = read_first(text, {.bonds = io::BondStrategy::templates});
        const auto& molecule = record.molecule;
        auto read_double_bonds = std::set<std::string>{};
        for (const auto& bond : molecule.bonds()) {
            if (bond.order() == chargefw::core::BondOrder::DOUBLE) {
                read_double_bonds.insert(
                    std::format("{}-{}", molecule.atom(bond.first_atom_index()).name(),
                                molecule.atom(bond.second_atom_index()).name()));
            }
        }
        CAPTURE(text);
        CHECK(read_double_bonds == double_bonds);
        for (std::size_t index = 0; index < molecule.atom_count(); ++index) {
            const auto& atom = molecule.atom(index);
            CHECK(atom.formal_charge() ==
                  (charges.contains(atom.name()) ? charges.at(atom.name()) : 0));
        }
    }
}

TEST_CASE("PDB templates connect nucleotide sugar atoms", "[io][pdb]") {
    // Deoxyadenosine monophosphate: the sugar joins the phosphate and the base, so all 21 heavy
    // atoms form one fragment with three rings.
    const auto record =
        read_first(pdb_residue("DA", "P OP1 OP2 O5' C5' C4' O4' C3' O3' C2' C1' N9 C8 N7 C5 C6 N6 "
                                     "N1 C2 N3 C4"),
                   {.bonds = io::BondStrategy::templates});
    CHECK(record.molecule.atom_count() == 21);
    CHECK(record.molecule.bond_count() == 23);
    CHECK(record.import_warnings.empty());
}

TEST_CASE("input options and source descriptions share stable string conversion", "[io][options]") {
    CHECK(io::record_selection_from_string("polymers-and-ligands") ==
          io::RecordSelection::polymers_and_ligands);
    CHECK(io::to_string(io::RecordSelection::polymers_and_ligands) == "polymers-and-ligands");
    CHECK(io::bond_strategy_from_string("explicit") == io::BondStrategy::explicit_bonds);
    CHECK(io::to_string(io::BondStrategy::hybrid) == "hybrid");
    CHECK(chargefw::io::conformer_selection_from_string("first") ==
          chargefw::io::ConformerSelection::first);
    CHECK(chargefw::io::to_string(chargefw::io::ConformerSelection::all) == "all");

    CHECK(chargefw::io::molecular_source_format_from_string("molecule-json") ==
          chargefw::io::MolecularSourceFormat::molecule_json);
    CHECK(chargefw::io::to_string(chargefw::io::MolecularSourceFormat::mmcif) == "mmcif");
    CHECK(chargefw::io::to_string(chargefw::io::SourceConnectivity::explicitly_empty) ==
          "explicitly-empty");

    CHECK_THROWS_AS(io::record_selection_from_string("unknown"), std::invalid_argument);
    CHECK_THROWS_AS(chargefw::io::molecular_source_format_from_string("cif"),
                    std::invalid_argument);
}

TEST_CASE("PDB input preserves selected model identity and connectivity", "[io][pdb]") {
    std::istringstream input{R"pdb(HEADER    TEST PDB
TITLE     TWO MODELS
MODEL        1
ATOM      1  O   HOH A   1       0.000   0.000   0.000  1.00 20.00           O  
ATOM      2  H1 AHOH A   1       0.957   0.000   0.000  1.00 20.00           H  
ATOM      3  H1 BHOH A   1       9.000   0.000   0.000  1.00 20.00           H  
ENDMDL
MODEL        2
ATOM      1  O   HOH A   1       0.100   0.000   0.000  1.00 20.00           O  
ATOM      2  H1 AHOH A   1       1.057   0.000   0.000  1.00 20.00           H  
ATOM      3  H1 BHOH A   1       9.100   0.000   0.000  1.00 20.00           H  
ENDMDL
CONECT    1    2
END
)pdb"};

    auto reader = io::PdbReader{input, "water.pdb"};
    const auto first = reader.next();
    REQUIRE(first.has_value());
    CHECK(first->identity.source_name == "water.pdb");
    CHECK(first->identity.record_index == 0);
    // Without a HEADER record, the molecule is named by the file.
    CHECK(first->identity.record_id.display_string() == "water");
    CHECK(first->molecule.atom_count() == 2);
    CHECK(first->molecule.bond_count() == 0);
    CHECK(first->molecule.atom(0).atomic_number() == 8);
    CHECK(first->molecule.atom(1).atomic_number() == 1);
    CHECK(first->molecule.atom(1).name() == "H1");
    CHECK(first->molecule.conformer_count() == 2);
    CHECK(first->molecule.conformer(0)[1].x == 0.957);

    CHECK(first->molecule.conformer(1).name() == "2");
    CHECK(first->molecule.conformer(1)[0].x == 0.1);
    CHECK(first->molecule.conformer(1)[1].x == 1.057);
    REQUIRE(first->import_metadata.has_value());
    const auto& mapping = *first->import_metadata;
    CHECK(mapping.format == chargefw::io::MolecularSourceFormat::pdb);
    CHECK(mapping.source_connectivity == chargefw::io::SourceConnectivity::present);
    CHECK(mapping.policy.alternate_locations == "first-source-order");
    REQUIRE(mapping.conformers.size() == 2);
    REQUIRE(mapping.components.size() == 1);
    CHECK(mapping.components[0] == chargefw::io::SourceComponentInstance{"HOH", {0, 1}});
    CHECK(mapping.conformers[0].id == "1");
    CHECK(mapping.conformers[1].id == "2");
    CHECK(mapping.conformers[0].sites[0].position == 0);
    CHECK(mapping.conformers[1].sites[0].position == 3);
    REQUIRE(mapping.atoms[1].structural_labels != nullptr);
    CHECK(mapping.atoms[1].structural_labels->author.atom == "H1");
    CHECK(mapping.atoms[1].structural_labels->author.residue == "HOH");
    CHECK(mapping.atoms[1].structural_labels->author.chain == "A");
    CHECK(mapping.atoms[1].structural_labels->alternate_location == "A");
    CHECK_FALSE(reader.next().has_value());

    {
        std::istringstream selection_input{
            R"pdb(ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  
HETATM    2  C1  LIG A   2       1.000   0.000   0.000  1.00 20.00           C  
HETATM    3  O   HOH A   3       2.000   0.000   0.000  1.00 20.00           O  
HETATM    4  O   WAT A   4       3.000   0.000   0.000  1.00 20.00           O
HETATM    5  O   HOH B   3       4.000   0.000   0.000  1.00 20.00           O
END
)pdb"};
        auto all_reader = io::PdbReader{selection_input};
        const auto all_record = all_reader.next();
        REQUIRE(all_record.has_value());
        CHECK(all_record->molecule.atom_count() == 5);
        REQUIRE(all_record->import_metadata.has_value());
        const auto& components = all_record->import_metadata->components;
        REQUIRE(components.size() == 5);
        CHECK(components[0] == chargefw::io::SourceComponentInstance{"ALA", {0}});
        CHECK(components[1] == chargefw::io::SourceComponentInstance{"LIG", {1}});
        CHECK(components[2] == chargefw::io::SourceComponentInstance{"HOH", {2}});
        CHECK(components[3] == chargefw::io::SourceComponentInstance{"WAT", {3}});
        CHECK(components[4] == chargefw::io::SourceComponentInstance{"HOH", {4}});
        CHECK(all_record->import_metadata->atoms[2].structural_labels->author.chain == "A");
        CHECK(all_record->import_metadata->atoms[4].structural_labels->author.chain == "B");

        const auto ligand_and_water_input =
            R"pdb(ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  
HETATM    2  C1  LIG A   2       1.000   0.000   0.000  1.00 20.00           C  
HETATM    3  O   HOH A   3       2.000   0.000   0.000  1.00 20.00           O  
HETATM    4  O   WAT A   4       3.000   0.000   0.000  1.00 20.00           O
END
)pdb";
        for (const auto selection :
             {io::RecordSelection::polymers_and_ligands, io::RecordSelection::polymers}) {
            CHECK(read_first(ligand_and_water_input, {.selection = selection})
                      .molecule.atom_count() == 2);
        }

        std::istringstream modified_polymer_input{
            R"pdb(ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C
HETATM    2 SE   MSE A   2       1.000   0.000   0.000  1.00 20.00          SE
TER
HETATM    3  C1  LIG A   3       2.000   0.000   0.000  1.00 20.00           C
END
)pdb"};
        auto modified_polymer_reader =
            io::PdbReader{modified_polymer_input, {}, {.selection = io::RecordSelection::polymers}};
        const auto modified_polymer_record = modified_polymer_reader.next();
        REQUIRE(modified_polymer_record.has_value());
        CHECK(modified_polymer_record->molecule.atom_count() == 2);
        CHECK(modified_polymer_record->molecule.atom(0).name() == "CA");
        CHECK(modified_polymer_record->molecule.atom(1).name() == "SE");
    }

    {
        const auto strategy_input =
            R"pdb(SSBOND   1 CYS A   3    CYS A   4                          
LINK         C   ALA A   1                 C1  LIG A   5
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N  
ATOM      2  CA  ALA A   1       1.450   0.000   0.000  1.00 20.00           C  
ATOM      3  C   ALA A   1       2.450   1.000   0.000  1.00 20.00           C  
ATOM      4  O   ALA A   1       3.450   1.000   0.000  1.00 20.00           O  
ATOM      5  CB  ALA A   1       1.450  -1.000   0.000  1.00 20.00           C  
ATOM      6  N   GLY A   2       2.200   2.200   0.000  1.00 20.00           N  
ATOM      7  CA  GLY A   2       3.200   3.200   0.000  1.00 20.00           C  
ATOM      8  C   GLY A   2       4.200   3.200   0.000  1.00 20.00           C  
ATOM      9  O   GLY A   2       5.200   3.200   0.000  1.00 20.00           O  
ATOM     10  SG  CYS A   3       6.200   3.200   0.000  1.00 20.00           S  
ATOM     11  SG  CYS A   4       7.200   3.200   0.000  1.00 20.00           S  
HETATM   12  C1  LIG A   5       8.200   3.200   0.000  1.00 20.00           C  
CONECT   11   12
END
)pdb";
        const auto read_strategy = [&](const io::BondStrategy strategy) {
            return read_first(strategy_input, {.bonds = strategy});
        };

        CHECK(read_strategy(io::BondStrategy::none).molecule.bond_count() == 0);
        const auto peptide = read_strategy(io::BondStrategy::templates);
        CHECK(peptide.molecule.bond_count() == 8);
        CHECK(read_strategy(io::BondStrategy::explicit_bonds).molecule.bond_count() == 2);
        CHECK(read_strategy(io::BondStrategy::hybrid).molecule.bond_count() == 10);

        REQUIRE(peptide.import_metadata.has_value());
        const auto& components = peptide.import_metadata->components;
        REQUIRE(components.size() == 5);
        CHECK(components[0].component_id == "ALA");
        CHECK(components[0].atom_indices.size() == 5);
        CHECK(components[1].component_id == "GLY");
        CHECK(components[1].atom_indices.size() == 4);
        bool peptide_link = false;
        for (const auto& bond : peptide.molecule.bonds()) {
            peptide_link = peptide_link ||
                           (bond.first_atom_index() == 2 && bond.second_atom_index() == 5) ||
                           (bond.first_atom_index() == 5 && bond.second_atom_index() == 2);
        }
        CHECK(peptide_link);
    }

    {
        const auto duplicate_input =
            R"pdb(ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N
ATOM      2  CA  ALA A   1       1.450   0.000   0.000  1.00 20.00           C
CONECT    1    2
END
)pdb";
        for (const auto strategy : {io::BondStrategy::templates, io::BondStrategy::explicit_bonds,
                                    io::BondStrategy::hybrid}) {
            CHECK(read_first(duplicate_input, {.bonds = strategy}).molecule.bond_count() == 1);
        }
    }

    {
        // The ligand pair comes first in CONECT but follows the template carbonyl in hybrid.
        const auto carbonyl_input =
            R"pdb(HETATM    1  C1  LIG A   1       5.000   0.000   0.000  1.00 20.00           C
HETATM    2  O1  LIG A   1       6.230   0.000   0.000  1.00 20.00           O
ATOM      3  C   ALA A   2       0.000   0.000   0.000  1.00 20.00           C
ATOM      4  O   ALA A   2       1.230   0.000   0.000  1.00 20.00           O
CONECT    1    2
CONECT    3    4
END
)pdb";
        const auto record = read_first(carbonyl_input, {.bonds = io::BondStrategy::hybrid});
        REQUIRE(record.molecule.bond_count() == 2);
        CHECK(record.molecule.bond(0).order() == chargefw::core::BondOrder::DOUBLE);
        REQUIRE(record.import_warnings.size() == 1);
        const auto& warning = record.import_warnings.front();
        CHECK(warning.kind == io::ImportIssueKind::bond_orders_unknown);
        REQUIRE(warning.location.bond_indices == std::vector<std::size_t>{1});
        const auto& ligand = record.molecule.bond(1);
        CHECK(ligand.first_atom_index() == 0);
        CHECK(ligand.second_atom_index() == 1);

        const auto explicit_record =
            read_first(carbonyl_input, {.bonds = io::BondStrategy::explicit_bonds});
        REQUIRE(explicit_record.molecule.bond_count() == 2);
        CHECK(explicit_record.molecule.bond(0).first_atom_index() == 0);
        CHECK(explicit_record.molecule.bond(1).order() == chargefw::core::BondOrder::SINGLE);
        REQUIRE(explicit_record.import_warnings.size() == 1);
        CHECK(explicit_record.import_warnings[0].kind == io::ImportIssueKind::bond_orders_unknown);
        CHECK(explicit_record.import_warnings[0].location.bond_indices ==
              std::vector<std::size_t>{0, 1});
    }
}

TEST_CASE("PDB input selects the first noncontiguous alternate location", "[io][pdb]") {
    std::istringstream input{
        R"pdb(HETATM    1  C1 BLIG A   1       0.000   0.000   0.000  1.00  0.00           C
HETATM    2  O1  LIG A   1       1.000   0.000   0.000  1.00  0.00           O
HETATM    3  C1 ALIG A   1       2.000   0.000   0.000  1.00  0.00           C
END
)pdb"};

    auto reader = io::PdbReader{input};
    const auto record = reader.next();
    REQUIRE(record.has_value());
    REQUIRE(record->molecule.atom_count() == 2);
    CHECK(record->molecule.atom(0).name() == "C1");
    CHECK(record->molecule.atom(1).name() == "O1");
    CHECK(record->molecule.conformer(0)[0].x == 0.0);
    CHECK(record->molecule.conformer(0)[1].x == 1.0);
    REQUIRE(record->import_metadata.has_value());
    REQUIRE(record->import_metadata->components.size() == 1);
    CHECK(record->import_metadata->components[0] ==
          chargefw::io::SourceComponentInstance{"LIG", {0, 1}});
    CHECK_FALSE(record->import_metadata->conformers[0].id.has_value());
    CHECK(record->import_metadata->atoms[0].position == 0);
    CHECK(record->import_metadata->atoms[0].id == "1");
    CHECK(record->import_metadata->atoms[1].position == 1);
    CHECK(record->import_metadata->atoms[1].id == "2");
    REQUIRE(record->import_metadata->atoms[0].structural_labels != nullptr);
    CHECK(record->import_metadata->atoms[0].structural_labels->alternate_location == "B");
}

TEST_CASE("PDB input is named by the ID code of its HEADER record", "[io][pdb]") {
    std::istringstream input{
        R"pdb(HEADER    HYDROLASE                               01-JAN-00   1ABC
HETATM    1  O   HOH A   1       0.000   0.000   0.000  1.00 20.00           O
END
)pdb"};
    auto reader = io::PdbReader{input, "0123abcd_upload.pdb"};
    const auto record = reader.next();
    REQUIRE(record.has_value());
    CHECK(record->identity.record_id.display_string() == "1ABC");
    CHECK(record->identity.source_name == "0123abcd_upload.pdb");
}

TEST_CASE("PDB input reads one residue of alternative residues at a position", "[io][pdb]") {
    // SER and PRO are alternatives at A 2 with coincident backbone atoms; GLY at 2A is a neighbor.
    const auto record = read_first(
        R"pdb(ATOM      1  N  ASER A   2       0.000   0.000   0.000  0.40 20.00           N
ATOM      2  CA ASER A   2       1.458   0.000   0.000  0.40 20.00           C
ATOM      3  CB ASER A   2       1.990   1.430   0.000  0.40 20.00           C
ATOM      4  N  CPRO A   2       0.000   0.000   0.000  0.60 20.00           N
ATOM      5  CA CPRO A   2       1.458   0.000   0.000  0.60 20.00           C
ATOM      6  CD CPRO A   2      -0.700   1.300   0.000  0.60 20.00           C
ATOM      7  N   GLY A   2A      3.000   0.000   0.000  1.00 20.00           N
END
)pdb",
        {.bonds = io::BondStrategy::templates});

    REQUIRE(record.import_metadata.has_value());
    const auto& metadata = *record.import_metadata;
    CHECK(metadata.components ==
          std::vector<io::SourceComponentInstance>{{"SER", {0, 1, 2}}, {"GLY", {3}}});
    const auto ids = metadata.atoms |
                     std::views::transform([](const auto& atom) { return atom.id.value_or(""); }) |
                     std::ranges::to<std::vector<std::string>>();
    CHECK(ids == std::vector<std::string>{"1", "2", "3", "7"});
    CHECK(metadata.atoms[0].structural_labels->alternate_location == "A");
    // N-CA and CA-CB of SER; the N-CD bond of PRO is omitted with its atoms.
    CHECK(record.molecule.bond_count() == 2);

    // Residues that repeat a number in another segment or without alternate locations are all
    // there; the second ligand has its atom at two locations.
    const auto repeated = read_first(
        R"pdb(HETATM    1  C1  LIG A   1       0.000   0.000   0.000  1.00 20.00      SEG1 C
HETATM    2  C1 ALIG A   1       5.000   0.000   0.000  0.50 20.00      SEG2 C
HETATM    3  C1 BLIG A   1       5.100   0.000   0.000  0.50 20.00      SEG2 C
HETATM    4  O   HOH A   1      10.000   0.000   0.000  1.00 20.00      SEG2 O
END
)pdb");
    const auto repeated_ids =
        repeated.import_metadata->atoms |
        std::views::transform([](const auto& atom) { return atom.id.value_or(""); }) |
        std::ranges::to<std::vector<std::string>>();
    CHECK(repeated_ids == std::vector<std::string>{"1", "2", "4"});
}

TEST_CASE("PDB input reports missing residue bonds and bond orders", "[io][pdb]") {
    using Strategy = io::BondStrategy;
    // A ligand without a template, and a single-atom ion that is unbonded by design.
    const std::string atoms =
        R"pdb(HETATM    1  C1  LIG A   1       0.000   0.000   0.000  1.00 20.00           C
HETATM    2  O1  LIG A   1       1.200   0.000   0.000  1.00 20.00           O
HETATM    3  C2  LIG A   1       0.000   1.400   0.000  1.00 20.00           C
HETATM    4  H1  LIG A   1      -1.000   0.000   0.000  1.00 20.00           H
HETATM    5 NA    NA A   2       5.000   0.000   0.000  1.00 20.00          NA
)pdb";
    const auto read_codes = [&atoms](const std::string& connectivity, const Strategy strategy) {
        auto codes = std::vector<std::string>{};
        for (const auto& issue :
             read_first(atoms + connectivity + "END\n", {.bonds = strategy}).import_warnings) {
            codes.emplace_back(chargefw::io::to_string(issue.kind));
        }
        return codes;
    };

    CHECK(read_codes("", Strategy::hybrid) == std::vector<std::string>{"unbonded_residue_atoms"});
    // The warning lists the ligand's atoms, not the ion.
    const auto unbonded = read_first(atoms + "END\n", {.bonds = Strategy::hybrid});
    REQUIRE(unbonded.import_warnings.size() == 1);
    CHECK(unbonded.import_warnings[0].location.atom_index == 0);
    CHECK(unbonded.import_warnings[0].location.atom_indices ==
          std::vector<std::size_t>{0, 1, 2, 3});
    CHECK(read_codes("", Strategy::none).empty());

    // Two bonds between ligand heavy atoms, one listed both ways, plus a bond to hydrogen, a link
    // to the ion, and a missing serial. The hydrogen bond is single for certain and not reported.
    const std::string connectivity = "CONECT    1    5    2    3    4   99\nCONECT    2    1\n";
    const auto bonded = read_first(atoms + connectivity + "END\n", {.bonds = Strategy::hybrid});
    REQUIRE(bonded.molecule.bond_count() == 4);
    REQUIRE(bonded.import_warnings.size() == 1);
    const auto& issue = bonded.import_warnings.front();
    CHECK(io::to_string(issue.kind) == "bond_orders_unknown");
    REQUIRE(issue.location.bond_indices.size() == 2);
    CHECK(issue.location.bond_indices[0] < issue.location.bond_indices[1]);
    for (const auto index : issue.location.bond_indices) {
        REQUIRE(index < bonded.molecule.bond_count());
        const auto& bond = bonded.molecule.bond(index);
        CHECK(bond.first_atom_index() < 3);
        CHECK(bond.second_atom_index() < 3);
        CHECK(bond.order() == chargefw::core::BondOrder::SINGLE);
    }
    CHECK(read_codes(connectivity, Strategy::none).empty());
    CHECK(read_codes(connectivity, Strategy::templates) ==
          std::vector<std::string>{"unbonded_residue_atoms"});
}

TEST_CASE("PDB explicit connections retain address wildcard behavior", "[io][pdb]") {
    const auto read_link = [](const std::string_view link) {
        auto record = read_first(std::string{link} + R"pdb(
ATOM      1  C  AALA A   1       0.000   0.000   0.000  1.00 20.00           C
HETATM    2  C1 ALIG A   5       1.000   0.000   0.000  1.00 20.00           C
END
)pdb",
                                 {.bonds = io::BondStrategy::explicit_bonds});
        CHECK(record.molecule.atom_count() == 2);
        return record;
    };

    CHECK(read_link("LINK         C   ALA A   1                 C1  LIG A   5")
              .molecule.bond_count() == 1);
    CHECK(read_link("LINK         C       A   1                 C1  LIG A   5")
              .molecule.bond_count() == 1);

    // A link to the ligand of a symmetry-related copy is not a bond of this molecule.
    const auto symmetry =
        read_link("LINK         C   ALA A   1                 C1  LIG A   5     1555   2555");
    CHECK(symmetry.molecule.bond_count() == 0);
    REQUIRE(symmetry.import_warnings.size() == 1);
    CHECK(symmetry.import_warnings.front().kind ==
          chargefw::io::ImportIssueKind::symmetry_connections_ignored);
}

TEST_CASE("PDB explicit connections require a unique residue for each bonded atom", "[io][pdb]") {
    // LINK has no segment ID. Only the second LIG has C1, but that does not identify the residue.
    const auto read_link = [](const std::string_view link) {
        return read_first(std::string{link} + R"pdb(
HETATM    1  O1  LIG A   1       0.000   0.000   0.000  1.00 20.00      SEG1 O
HETATM    2  C1 ALIG A   1       5.000   0.000   0.000  1.00 20.00      SEG2 C
ATOM      3  C   ALA A   5       6.000   0.000   0.000  1.00 20.00           C
END
)pdb",
                          {.bonds = io::BondStrategy::explicit_bonds});
    };
    for (const auto* link : {"LINK         C   ALA A   9                 C1  LIG A   1",
                             "LINK         C1  LIG A   1                 C   ALA A   9",
                             "LINK         C   ALA A   5                 C1 BLIG A   1"}) {
        CAPTURE(link);
        CHECK(read_link(link).molecule.bond_count() == 0);
    }
    const auto symmetry =
        read_link("LINK         C   ALA A   5                 C1  LIG A   1     1555   2555");
    CHECK(symmetry.molecule.bond_count() == 0);
    REQUIRE(symmetry.import_warnings.size() == 1);
    CHECK(symmetry.import_warnings.front().kind ==
          io::ImportIssueKind::symmetry_connections_ignored);

    const auto error = chargefw::test::thrown<io::ParseError>([&] {
        static_cast<void>(read_link("LINK         C   ALA A   5                 C1  LIG A   1"));
    });
    REQUIRE(error.has_value());
    CHECK(error->kind() == io::ParseErrorKind::malformed);
    CHECK(error->reason().contains("LIG A 1"));
    CHECK(error->reason().contains("C1"));
}

TEST_CASE("PDB source mapping follows Gemmi record handling", "[io][pdb]") {
    std::istringstream input{R"pdb(mOdEl        1
aToM      1  C1  LIG A   1       0.000   0.000   0.000  1.00 20.00           C
eNdMdL
eNd
cOnEcT    1    2
aToM      2  O1  LIG A   1       1.000   0.000   0.000  1.00 20.00           O
)pdb"};

    auto reader = io::PdbReader{input};
    const auto record = reader.next();
    REQUIRE(record.has_value());
    CHECK(record->molecule.atom_count() == 1);
    REQUIRE(record->import_metadata.has_value());
    const auto& mapping = *record->import_metadata;
    CHECK(mapping.source_connectivity == chargefw::io::SourceConnectivity::absent);
    REQUIRE(mapping.conformers.size() == 1);
    CHECK(mapping.conformers[0].id == "1");
    REQUIRE(mapping.conformers[0].sites.size() == 1);
    CHECK(mapping.conformers[0].sites[0].position == 0);
}

TEST_CASE("PDB mapping retains model-specific alternate locations", "[io][pdb]") {
    std::istringstream input{R"pdb(MODEL        1
ATOM      1  C1 ALIG A   1       0.000   0.000   0.000  1.00 20.00           C
ENDMDL
MODEL        2
ATOM      1  C1 BLIG A   1       0.100   0.000   0.000  1.00 20.00           C
ENDMDL
END
)pdb"};
    auto reader = io::PdbReader{input};
    const auto record = reader.next();
    REQUIRE(record.has_value());
    REQUIRE(record->import_metadata.has_value());
    const auto& conformers = record->import_metadata->conformers;
    REQUIRE(conformers.size() == 2);
    REQUIRE(conformers[0].sites[0].structural_labels != nullptr);
    REQUIRE(conformers[1].sites[0].structural_labels != nullptr);
    CHECK(conformers[0].sites[0].structural_labels->alternate_location == "A");
    CHECK(conformers[1].sites[0].structural_labels->alternate_location == "B");
}

TEST_CASE("PDB input rejects an atom repeated at one alternate location", "[io][pdb]") {
    // Without alternate locations, and at B after A, which is the one read.
    for (const std::string_view locations : {"  ", "ABB"}) {
        auto text = std::string{};
        for (std::size_t index = 0; index < locations.size(); ++index) {
            text +=
                std::format("ATOM  {:>5}  C1 {}LIG A   1       {}.000   0.000   0.000  1.00 20.00"
                            "           C\n",
                            index + 1, locations[index], index);
        }
        CAPTURE(text);
        std::istringstream input{text + "END\n"};
        const auto error = chargefw::test::thrown<io::ParseError>(
            [&] { static_cast<void>(io::PdbReader{input}); });
        REQUIRE(error.has_value());
        CHECK(error->kind() == io::ParseErrorKind::malformed);
        CHECK(error->reason().contains("LIG A 1"));
        CHECK(error->reason().contains("C1"));
        if (locations.back() != ' ') {
            CHECK(error->reason().ends_with(locations.back()));
        }
    }
}

TEST_CASE("PDB mapping retains source labels through Gemmi normalization", "[io][pdb]") {
    std::istringstream input{
        R"pdb(ATOM      1  C1  LIGAB0001      0.000   0.000   0.000  1.00 20.00           C
ATOM      2  O1  LIGABA000      1.000   0.000   0.000  1.00 20.00           O
END
)pdb"};
    auto reader = io::PdbReader{input};
    const auto record = reader.next();
    REQUIRE(record.has_value());
    REQUIRE(record->import_metadata.has_value());
    const auto& atoms = record->import_metadata->atoms;
    REQUIRE(atoms.size() == 2);
    REQUIRE(atoms[0].structural_labels != nullptr);
    REQUIRE(atoms[1].structural_labels != nullptr);
    CHECK(atoms[0].structural_labels->author.chain == "AB");
    CHECK(atoms[0].structural_labels->author.sequence == "0001");
    CHECK(atoms[1].structural_labels->author.chain == "AB");
    CHECK(atoms[1].structural_labels->author.sequence == "A000");
}

TEST_CASE("PDB input rejects empty and incompatible selected models", "[io][pdb]") {
    for (const auto* const text : {"", "HEADER    TEST\nEND\n", "MODEL        1\nENDMDL\nEND\n"}) {
        CAPTURE(text);
        std::istringstream input{text};
        const auto error = chargefw::test::thrown<io::ParseError>(
            [&] { static_cast<void>(io::PdbReader{input, "input.pdb", {}}); });
        REQUIRE(error.has_value());
        CHECK(error->location().source_name == "input.pdb");
    }

    {
        std::istringstream input{R"pdb(MODEL        1
ATOM      1  C   UNL A   1       0.000   0.000   0.000  1.00 20.00           C
ENDMDL
MODEL        2
ATOM      1  O   UNL A   1       0.000   0.000   0.000  1.00 20.00           O
ENDMDL
END
)pdb"};
        CHECK_THROWS_AS(io::PdbReader{input}, std::exception);
    }

    {
        std::istringstream input{R"pdb(MODEL        1
ATOM      1  C1  UNL A   1       0.000   0.000   0.000  1.00 20.00           C
ENDMDL
MODEL        2
ATOM      1  C1  UNL B   1       0.100   0.000   0.000  1.00 20.00           C
ENDMDL
END
)pdb"};
        CHECK_THROWS_AS(io::PdbReader{input}, std::runtime_error);
    }
}

TEST_CASE("PDB template links respect chains and insertion codes", "[io][pdb]") {
    const auto read_bond_count = [](const std::string& text) {
        return read_first(text, {.bonds = io::BondStrategy::templates}).molecule.bond_count();
    };

    CHECK(read_bond_count(
              R"pdb(ATOM      1  C   ALA A  82       0.000   0.000   0.000  1.00 20.00           C
ATOM      2  N   GLY A  82A      1.000   0.000   0.000  1.00 20.00           N
ATOM      3  C   GLY A  82A      2.000   0.000   0.000  1.00 20.00           C
ATOM      4  N   SER A  82B      3.000   0.000   0.000  1.00 20.00           N
END
)pdb") == 2);

    CHECK(read_bond_count(
              R"pdb(ATOM      1  C   ALA A  82A      0.000   0.000   0.000  1.00 20.00           C
ATOM      2  N   GLY B  82B      1.000   0.000   0.000  1.00 20.00           N
END
)pdb") == 0);

    CHECK(read_bond_count(
              R"pdb(ATOM      1  C   ALA A  81       0.000   0.000   0.000  1.00 20.00           C
ATOM      2  N   GLY A  83       1.000   0.000   0.000  1.00 20.00           N
END
)pdb") == 0);

    CHECK(read_bond_count(
              R"pdb(ATOM      1  C   ALA A  82A      0.000   0.000   0.000  1.00 20.00           C
ATOM      2  N   GLY A  82C      1.000   0.000   0.000  1.00 20.00           N
END
)pdb") == 0);
}

TEST_CASE("PDB input gives the kind of refusal and passes a failed read on", "[io][pdb]") {
    const auto kind = [](const std::string& text, const io::StructureOptions options = {}) {
        std::istringstream input{text};
        const auto error = chargefw::test::thrown<io::ParseError>(
            [&] { static_cast<void>(io::PdbReader{input, "input.pdb", options}); });
        REQUIRE(error.has_value());
        CHECK(error->location().source_name == "input.pdb");
        return error->kind();
    };
    const auto atom = [](const std::string_view name, const std::string_view element,
                         const std::string_view record = "ATOM  ") {
        return std::format("{}    1 {:<4} HOH A   1       0.000   0.000   0.000  1.00  0.00"
                           "          {:>2}\n",
                           record, name, element);
    };
    CHECK(kind("") == io::ParseErrorKind::no_records);
    CHECK(kind(atom(" O", "O", "HETATM"), {.selection = io::RecordSelection::polymers}) ==
          io::ParseErrorKind::no_records);
    CHECK(kind("REMARK   1 Caf\xe9\n" + atom(" O", "O")) == io::ParseErrorKind::invalid_encoding);
    CHECK(kind(atom(" X", "Xx")) == io::ParseErrorKind::malformed);
    // Refused by the parsing library, and atoms that cannot be told apart.
    CHECK(kind("MODEL        1\n" + atom(" O", "O") + "MODEL        2\n" + atom(" O", "O") +
               "ENDMDL\nEND\n") == io::ParseErrorKind::malformed);
    CHECK(kind(atom(" O", "O") + atom(" O", "O")) == io::ParseErrorKind::malformed);
    CHECK(kind("MODEL        1\n" + atom(" C", "C") + "ENDMDL\nMODEL        2\n" + atom(" O", "O") +
               "ENDMDL\nEND\n") == io::ParseErrorKind::inconsistent_models);

    for (const auto state : {std::ios::failbit, std::ios::badbit}) {
        std::istringstream input{atom(" O", "O")};
        input.setstate(state);
        auto refused = false;
        auto failed = false;
        try {
            static_cast<void>(io::PdbReader{input, "input.pdb", {}});
        } catch (const io::ParseError&) {
            refused = true;
        } catch (const std::runtime_error&) {
            failed = true;
        }
        CHECK(failed);
        CHECK_FALSE(refused);
    }
}

TEST_CASE("PDB input rejects non-ASCII text", "[io][pdb]") {
    const auto read = [](const std::string& text) {
        std::istringstream input{text};
        static_cast<void>(io::PdbReader{input, {}, {}});
    };
    const auto error = chargefw::test::thrown<io::ParseError>(
        [&] { read("HEADER    TEST\nREMARK   1 Caf\xe9\n" + pdb_residue("HOH", "O")); });
    REQUIRE(error.has_value());
    const auto* location = &error->location();
    CHECK(location->line == 2);
    CHECK_FALSE(location->record_index.has_value());
}

TEST_CASE("A component dictionary gives PDB residues the bonds that nothing else does",
          "[io][pdb][components]") {
    // The ligand entry lacks N1. The alanine entry differs from the built-in template: a single
    // carbonyl bond, a hydrogen HZ9 on CB that the template does not have, and no CA.
    std::istringstream dictionary_text{"chargefw-ccd 1 2026.11.13\n"
                                       "LIG 3 2\n"
                                       "C1 C 0\n"
                                       "O1 O 0\n"
                                       "C2 C 0\n"
                                       "0 1 2\n"
                                       "0 2 1\n"
                                       "ALA 4 2\n"
                                       "C C 0\n"
                                       "O O 0\n"
                                       "CB C 0\n"
                                       "HZ9 H 0\n"
                                       "0 1 1\n"
                                       "2 3 1\n"};
    const auto dictionary = std::make_shared<const io::ComponentDictionary>(dictionary_text);
    const auto input =
        R"pdb(HETATM    1  C1  LIG A   1       5.000   0.000   0.000  1.00 20.00           C
HETATM    2  O1  LIG A   1       6.230   0.000   0.000  1.00 20.00           O
HETATM    3  C2  LIG A   1       5.000   1.500   0.000  1.00 20.00           C
HETATM    4  N1  LIG A   1       5.000   3.000   0.000  1.00 20.00           N
ATOM      5  C   ALA A   2       0.000   0.000   0.000  1.00 20.00           C
ATOM      6  O   ALA A   2       1.230   0.000   0.000  1.00 20.00           O
ATOM      7  CB  ALA A   2       0.000   3.000   0.000  1.00 20.00           C
ATOM      8  HZ9 ALA A   2       0.000   4.000   0.000  1.00 20.00           H
ATOM      9  CA  ALA A   2       0.000   1.500   0.000  1.00 20.00           C
CONECT    1    2
CONECT    3    4
END
)pdb";
    using chargefw::core::BondOrder;
    // The order of the bond between two atoms, or nothing when they are not bonded.
    const auto order = [](const io::ImportedMoleculeRecord& record, const std::size_t first,
                          const std::size_t second) -> std::optional<BondOrder> {
        for (const auto& bond : record.molecule.bonds()) {
            if ((bond.first_atom_index() == first && bond.second_atom_index() == second) ||
                (bond.first_atom_index() == second && bond.second_atom_index() == first)) {
                return bond.order();
            }
        }
        return std::nullopt;
    };
    const auto unknown_orders = [](const io::ImportedMoleculeRecord& record) {
        for (const auto& warning : record.import_warnings) {
            if (warning.kind == io::ImportIssueKind::bond_orders_unknown) {
                return warning.location.bond_indices;
            }
        }
        return std::vector<std::size_t>{};
    };
    // The atoms reported as having no bonds.
    const auto unbonded_atoms = [](const io::ImportedMoleculeRecord& record) {
        for (const auto& warning : record.import_warnings) {
            if (warning.kind == io::ImportIssueKind::unbonded_residue_atoms) {
                return warning.location.atom_indices;
            }
        }
        return std::vector<std::size_t>{};
    };

    {
        // Hybrid fills gaps bond by bond. The dictionary gives the order of one CONECT pair and
        // adds C1-C2; the other CONECT pair stays the only one of unknown order.
        const auto record =
            read_first(input, {.bonds = io::BondStrategy::hybrid, .components = dictionary});
        CHECK(order(record, 0, 1) == BondOrder::DOUBLE);
        CHECK(order(record, 0, 2) == BondOrder::SINGLE);
        CHECK(order(record, 2, 3) == BondOrder::SINGLE);
        const auto unknown = unknown_orders(record);
        REQUIRE(unknown.size() == 1);
        const auto& bond = record.molecule.bond(unknown.front());
        CHECK(bond.first_atom_index() == 2);
        CHECK(bond.second_atom_index() == 3);
        // The built-in template keeps its carbonyl order, and the dictionary adds the pair the
        // template lacks.
        CHECK(order(record, 4, 5) == BondOrder::DOUBLE);
        CHECK(order(record, 6, 7) == BondOrder::SINGLE);
        // Every atom has a bond, so nothing is reported: N1 and the alanine CA are not in their
        // entries, but the file bonds one and the built-in template the other.
        CHECK(unbonded_atoms(record).empty());
        REQUIRE(record.import_metadata.has_value());
        CHECK(record.import_metadata->policy.components == "2026.11.13");
    }
    {
        // A bond the file states with an order keeps it: a LINK record is a single bond, although
        // the dictionary has the pair as double.
        const auto linked =
            std::string{"LINK         C1  LIG A   1                 O1  LIG A   1\n"} + input;
        const auto record =
            read_first(linked, {.bonds = io::BondStrategy::hybrid, .components = dictionary});
        CHECK(order(record, 0, 1) == BondOrder::SINGLE);
        CHECK(order(record, 0, 2) == BondOrder::SINGLE);
        CHECK(unknown_orders(record).size() == 1);
    }
    {
        // Templates take one template per residue: the dictionary's for the ligand, the built-in
        // one alone for alanine. CONECT is not read.
        const auto record =
            read_first(input, {.bonds = io::BondStrategy::templates, .components = dictionary});
        CHECK(order(record, 0, 1) == BondOrder::DOUBLE);
        CHECK(order(record, 0, 2) == BondOrder::SINGLE);
        CHECK_FALSE(order(record, 2, 3).has_value());
        CHECK(order(record, 4, 5) == BondOrder::DOUBLE);
        CHECK_FALSE(order(record, 6, 7).has_value());
        // The atoms that no template names are left without bonds and reported.
        CHECK(unbonded_atoms(record) == std::vector<std::size_t>{3, 7});
    }
    {
        // Explicit bonds are those of the file, whatever the dictionary has.
        const auto record = read_first(
            input, {.bonds = io::BondStrategy::explicit_bonds, .components = dictionary});
        CHECK(record.molecule.bond_count() == 2);
        CHECK(order(record, 0, 1) == BondOrder::SINGLE);
        CHECK(unknown_orders(record).size() == 2);
        CHECK_FALSE(record.import_metadata->policy.components.has_value());
    }
    {
        // The record names the dictionary it was read with, also when the dictionary adds nothing:
        // the alanine carbonyl is the template's.
        const auto carbonyl =
            R"pdb(ATOM      1  C   ALA A   2       0.000   0.000   0.000  1.00 20.00           C
ATOM      2  O   ALA A   2       1.230   0.000   0.000  1.00 20.00           O
END
)pdb";
        const auto record =
            read_first(carbonyl, {.bonds = io::BondStrategy::hybrid, .components = dictionary});
        CHECK(order(record, 0, 1) == BondOrder::DOUBLE);
        CHECK(record.import_metadata->policy.components == "2026.11.13");
    }
    {
        // Without a dictionary, hybrid reads the file as before.
        const auto record = read_first(input, {.bonds = io::BondStrategy::hybrid});
        CHECK(order(record, 0, 1) == BondOrder::SINGLE);
        CHECK_FALSE(order(record, 0, 2).has_value());
        CHECK_FALSE(order(record, 6, 7).has_value());
        CHECK(unknown_orders(record).size() == 2);
    }
}
