#include <chargefw/adapters/conformer_selection.h>
#include <chargefw/adapters/gemmi/input_options.h>
#include <chargefw/adapters/gemmi/pdb_input.h>
#include <chargefw/core/bond.h>
#include <snitch/snitch.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <exception>
#include <format>
#include <map>
#include <ranges>
#include <set>
#include <sstream>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <vector>

namespace gemmi_adapter = chargefw::adapters::gemmi;
namespace pdb = gemmi_adapter::pdb_input;

namespace {

[[nodiscard]] auto read_first(const std::string& text,
                              const gemmi_adapter::InputOptions& options = {}) {
    std::istringstream input{text};
    auto reader = pdb::PdbReader{input, {}, options};
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

static_assert(!std::is_copy_constructible_v<pdb::PdbReader> &&
              !std::is_copy_assignable_v<pdb::PdbReader>);
static_assert(std::is_move_constructible_v<pdb::PdbReader>);

TEST_CASE("PDB templates supplement N-terminal hydrogen connectivity", "[adapters][pdb]") {
    using Strategy = gemmi_adapter::BondStrategy;
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
        auto reader = pdb::PdbReader{input, {}, {.bond_strategy = strategy}};
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
        {.bond_strategy = Strategy::templates});
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
        {.bond_strategy = Strategy::templates});
    CHECK(std::ranges::any_of(unlinked.molecule.bonds(), [](const auto& bond) {
        return bond.first_atom_index() == 2 && bond.second_atom_index() == 4;
    }));
}

TEST_CASE("PDB templates follow the protonation state of histidine and arginine",
          "[adapters][pdb]") {
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
        const auto record =
            read_first(text, {.bond_strategy = gemmi_adapter::BondStrategy::templates});
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

TEST_CASE("PDB templates connect nucleotide sugar atoms", "[adapters][pdb]") {
    // Deoxyadenosine monophosphate: the sugar joins the phosphate and the base, so all 21 heavy
    // atoms form one fragment with three rings.
    const auto record =
        read_first(pdb_residue("DA", "P OP1 OP2 O5' C5' C4' O4' C3' O3' C2' C1' N9 C8 N7 C5 C6 N6 "
                                     "N1 C2 N3 C4"),
                   {.bond_strategy = gemmi_adapter::BondStrategy::templates});
    CHECK(record.molecule.atom_count() == 21);
    CHECK(record.molecule.bond_count() == 23);
    CHECK(record.issues.empty());
}

TEST_CASE("structural input options share stable string conversion", "[adapters][options]") {
    CHECK(gemmi_adapter::record_selection_from_string("polymers-and-ligands") ==
          gemmi_adapter::RecordSelection::polymers_and_ligands);
    CHECK(gemmi_adapter::to_string(gemmi_adapter::RecordSelection::polymers_and_ligands) ==
          "polymers-and-ligands");
    CHECK(gemmi_adapter::bond_strategy_from_string("explicit") ==
          gemmi_adapter::BondStrategy::explicit_bonds);
    CHECK(gemmi_adapter::to_string(gemmi_adapter::BondStrategy::hybrid) == "hybrid");
    CHECK(chargefw::adapters::conformer_selection_from_string("first") ==
          chargefw::adapters::ConformerSelection::first);
    CHECK(chargefw::adapters::to_string(chargefw::adapters::ConformerSelection::all) == "all");

    CHECK_THROWS_AS(gemmi_adapter::record_selection_from_string("unknown"), std::invalid_argument);
}

TEST_CASE("PDB input preserves selected model identity and connectivity", "[adapters][pdb]") {
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

    auto reader = pdb::PdbReader{input, "water.pdb"};
    const auto first = reader.next();
    REQUIRE(first.has_value());
    CHECK(first->identity.source == "water.pdb");
    CHECK(first->identity.record_index == 0);
    CHECK_FALSE(first->identity.record_id.empty());
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
    CHECK(mapping.format == chargefw::adapters::MolecularSourceFormat::pdb);
    CHECK(mapping.source_connectivity == chargefw::adapters::SourceConnectivity::present);
    CHECK(mapping.alternate_location_selection == "first-source-order");
    REQUIRE(mapping.conformers.size() == 2);
    REQUIRE(mapping.components.size() == 1);
    CHECK(mapping.components[0] == chargefw::adapters::SourceComponentInstance{"HOH", {0, 1}});
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
        auto all_reader = pdb::PdbReader{selection_input};
        const auto all_record = all_reader.next();
        REQUIRE(all_record.has_value());
        CHECK(all_record->molecule.atom_count() == 5);
        REQUIRE(all_record->import_metadata.has_value());
        const auto& components = all_record->import_metadata->components;
        REQUIRE(components.size() == 5);
        CHECK(components[0] == chargefw::adapters::SourceComponentInstance{"ALA", {0}});
        CHECK(components[1] == chargefw::adapters::SourceComponentInstance{"LIG", {1}});
        CHECK(components[2] == chargefw::adapters::SourceComponentInstance{"HOH", {2}});
        CHECK(components[3] == chargefw::adapters::SourceComponentInstance{"WAT", {3}});
        CHECK(components[4] == chargefw::adapters::SourceComponentInstance{"HOH", {4}});
        CHECK(all_record->import_metadata->atoms[2].structural_labels->author.chain == "A");
        CHECK(all_record->import_metadata->atoms[4].structural_labels->author.chain == "B");

        const auto ligand_and_water_input =
            R"pdb(ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  
HETATM    2  C1  LIG A   2       1.000   0.000   0.000  1.00 20.00           C  
HETATM    3  O   HOH A   3       2.000   0.000   0.000  1.00 20.00           O  
HETATM    4  O   WAT A   4       3.000   0.000   0.000  1.00 20.00           O
END
)pdb";
        for (const auto selection : {gemmi_adapter::RecordSelection::polymers_and_ligands,
                                     gemmi_adapter::RecordSelection::polymers}) {
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
        auto modified_polymer_reader = pdb::PdbReader{
            modified_polymer_input, {}, {.selection = gemmi_adapter::RecordSelection::polymers}};
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
        const auto read_strategy = [&](const gemmi_adapter::BondStrategy strategy) {
            return read_first(strategy_input, {.bond_strategy = strategy});
        };

        CHECK(read_strategy(gemmi_adapter::BondStrategy::none).molecule.bond_count() == 0);
        const auto peptide = read_strategy(gemmi_adapter::BondStrategy::templates);
        CHECK(peptide.molecule.bond_count() == 8);
        CHECK(read_strategy(gemmi_adapter::BondStrategy::explicit_bonds).molecule.bond_count() ==
              2);
        CHECK(read_strategy(gemmi_adapter::BondStrategy::hybrid).molecule.bond_count() == 10);

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
        for (const auto strategy :
             {gemmi_adapter::BondStrategy::templates, gemmi_adapter::BondStrategy::explicit_bonds,
              gemmi_adapter::BondStrategy::hybrid}) {
            CHECK(read_first(duplicate_input, {.bond_strategy = strategy}).molecule.bond_count() ==
                  1);
        }
    }

    {
        // CONECT carries no order, so hybrid keeps the template's carbonyl double bond.
        const auto carbonyl_input =
            R"pdb(ATOM      1  C   ALA A   1       0.000   0.000   0.000  1.00 20.00           C
ATOM      2  O   ALA A   1       1.230   0.000   0.000  1.00 20.00           O
CONECT    1    2
END
)pdb";
        const auto record =
            read_first(carbonyl_input, {.bond_strategy = gemmi_adapter::BondStrategy::hybrid});
        REQUIRE(record.molecule.bond_count() == 1);
        CHECK(record.molecule.bond(0).order() == chargefw::core::BondOrder::DOUBLE);
    }
}

TEST_CASE("PDB input selects the first noncontiguous alternate location", "[adapters][pdb]") {
    std::istringstream input{
        R"pdb(HETATM    1  C1 BLIG A   1       0.000   0.000   0.000  1.00  0.00           C
HETATM    2  O1  LIG A   1       1.000   0.000   0.000  1.00  0.00           O
HETATM    3  C1 ALIG A   1       2.000   0.000   0.000  1.00  0.00           C
END
)pdb"};

    auto reader = pdb::PdbReader{input};
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
          chargefw::adapters::SourceComponentInstance{"LIG", {0, 1}});
    CHECK_FALSE(record->import_metadata->conformers[0].id.has_value());
    CHECK(record->import_metadata->atoms[0].position == 0);
    CHECK(record->import_metadata->atoms[0].id == "1");
    CHECK(record->import_metadata->atoms[1].position == 1);
    CHECK(record->import_metadata->atoms[1].id == "2");
    REQUIRE(record->import_metadata->atoms[0].structural_labels != nullptr);
    CHECK(record->import_metadata->atoms[0].structural_labels->alternate_location == "B");
}

TEST_CASE("PDB input reports residue atoms without bonds", "[adapters][pdb]") {
    using Strategy = gemmi_adapter::BondStrategy;
    // A ligand without a template, and a single-atom ion that is unbonded by design.
    const std::string atoms =
        R"pdb(HETATM    1  C1  LIG A   1       0.000   0.000   0.000  1.00 20.00           C
HETATM    2  O1  LIG A   1       1.200   0.000   0.000  1.00 20.00           O
HETATM    3 NA    NA A   2       5.000   0.000   0.000  1.00 20.00          NA
)pdb";
    const auto read_codes = [&atoms](const std::string& connectivity, const Strategy strategy) {
        auto codes = std::vector<std::string>{};
        for (const auto& issue :
             read_first(atoms + connectivity + "END\n", {.bond_strategy = strategy}).issues) {
            codes.emplace_back(chargefw::adapters::to_string(issue.kind));
        }
        return codes;
    };

    CHECK(read_codes("", Strategy::hybrid) == std::vector<std::string>{"unbonded_residue_atoms"});
    // The warning lists the ligand's atoms, not the ion.
    const auto unbonded = read_first(atoms + "END\n", {.bond_strategy = Strategy::hybrid});
    REQUIRE(unbonded.issues.size() == 1);
    CHECK(unbonded.issues[0].location.atom_index == 0);
    CHECK(unbonded.issues[0].location.atom_indices == std::vector<std::size_t>{0, 1});
    CHECK(read_codes("CONECT    1    2\n", Strategy::hybrid).empty());
    CHECK(read_codes("", Strategy::none).empty());
}

TEST_CASE("PDB explicit connections retain address wildcard behavior", "[adapters][pdb]") {
    const auto read_link = [](const std::string_view link) {
        auto record = read_first(std::string{link} + R"pdb(
ATOM      1  C  AALA A   1       0.000   0.000   0.000  1.00 20.00           C
HETATM    2  C1 ALIG A   5       1.000   0.000   0.000  1.00 20.00           C
END
)pdb",
                                 {.bond_strategy = gemmi_adapter::BondStrategy::explicit_bonds});
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
    REQUIRE(symmetry.issues.size() == 1);
    CHECK(symmetry.issues.front().kind ==
          chargefw::adapters::ImportIssueKind::symmetry_connections_ignored);
}

TEST_CASE("PDB source mapping follows Gemmi record handling", "[adapters][pdb]") {
    std::istringstream input{R"pdb(mOdEl        1
aToM      1  C1  LIG A   1       0.000   0.000   0.000  1.00 20.00           C
eNdMdL
eNd
cOnEcT    1    2
aToM      2  O1  LIG A   1       1.000   0.000   0.000  1.00 20.00           O
)pdb"};

    auto reader = pdb::PdbReader{input};
    const auto record = reader.next();
    REQUIRE(record.has_value());
    CHECK(record->molecule.atom_count() == 1);
    REQUIRE(record->import_metadata.has_value());
    const auto& mapping = *record->import_metadata;
    CHECK(mapping.source_connectivity == chargefw::adapters::SourceConnectivity::absent);
    REQUIRE(mapping.conformers.size() == 1);
    CHECK(mapping.conformers[0].id == "1");
    REQUIRE(mapping.conformers[0].sites.size() == 1);
    CHECK(mapping.conformers[0].sites[0].position == 0);
}

TEST_CASE("PDB mapping retains model-specific alternate locations", "[adapters][pdb]") {
    std::istringstream input{R"pdb(MODEL        1
ATOM      1  C1 ALIG A   1       0.000   0.000   0.000  1.00 20.00           C
ENDMDL
MODEL        2
ATOM      1  C1 BLIG A   1       0.100   0.000   0.000  1.00 20.00           C
ENDMDL
END
)pdb"};
    auto reader = pdb::PdbReader{input};
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

TEST_CASE("PDB input rejects ambiguous source identities", "[adapters][pdb]") {
    std::istringstream input{
        R"pdb(ATOM      1  C1  LIG A   1       0.000   0.000   0.000  1.00 20.00           C
ATOM      2  C1  LIG A   1       0.100   0.000   0.000  1.00 20.00           C
END
)pdb"};

    CHECK_THROWS_AS(pdb::PdbReader{input}, std::runtime_error);
}

TEST_CASE("PDB mapping retains source labels through Gemmi normalization", "[adapters][pdb]") {
    std::istringstream input{
        R"pdb(ATOM      1  C1  LIGAB0001      0.000   0.000   0.000  1.00 20.00           C
ATOM      2  O1  LIGABA000      1.000   0.000   0.000  1.00 20.00           O
END
)pdb"};
    auto reader = pdb::PdbReader{input};
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

TEST_CASE("PDB input rejects empty and incompatible selected models", "[adapters][pdb]") {
    {
        std::istringstream input{"END\n"};
        CHECK_THROWS_AS(pdb::PdbReader{input}, std::exception);
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
        CHECK_THROWS_AS(pdb::PdbReader{input}, std::exception);
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
        CHECK_THROWS_AS(pdb::PdbReader{input}, std::runtime_error);
    }
}

TEST_CASE("PDB template links respect chains and insertion codes", "[adapters][pdb]") {
    const auto read_bond_count = [](const std::string& text) {
        return read_first(text, {.bond_strategy = gemmi_adapter::BondStrategy::templates})
            .molecule.bond_count();
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
