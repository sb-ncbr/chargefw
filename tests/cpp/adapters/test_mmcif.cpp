#include <chargefw/adapters/gemmi/mmcif_input.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/molecule.h>
#include <snitch/snitch.hpp>

#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>

namespace mmcif = chargefw::adapters::gemmi::mmcif_input;
namespace gemmi_adapter = chargefw::adapters::gemmi;

namespace {

// mmCIF data block opening with the atom-site loop header shared by these fixtures; each fixture
// appends its own atom rows.
[[nodiscard]] auto atom_site_header(const std::string_view block_name) -> std::string {
    return "data_" + std::string{block_name} + "\nloop_\n" + R"cif(_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num)cif";
}

[[nodiscard]] auto read_first(const std::string& text, const gemmi_adapter::InputOptions& options)
    -> chargefw::core::Molecule {
    std::istringstream input{text};
    auto reader = mmcif::MmcifReader{input, {}, options};
    auto record = reader.next();
    REQUIRE(record.has_value());
    return std::move(record->molecule);
}

} // namespace

static_assert(!std::is_copy_constructible_v<mmcif::MmcifReader> &&
              !std::is_copy_assignable_v<mmcif::MmcifReader>);
static_assert(std::is_move_constructible_v<mmcif::MmcifReader>);

TEST_CASE("mmCIF input preserves arbitrary atom-site IDs and label namespaces",
          "[adapters][mmcif]") {
    const auto header = R"cif(data_ids
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.label_entity_id
_atom_site.pdbx_PDB_model_num
)cif";
    std::istringstream input{
        std::string{header} +
        "HETATM Csite C L1 . LIG LA 7 ? 0 0 0 1 20 0 17 AUTH AC A1 E1 1\n"
        "HETATM 001 O L2 B LIG LA 7 ? 1 0 0 1 20 0 17 AUTH AC A2 E1 1\n"
        "HETATM 9223372036854775808 N L3 . LIG LA 7 ? 2 0 0 1 20 0 17 AUTH AC A3 E1 1\n"
        "HETATM '.' C L4 . LIG LA 7 ? 3 0 0 1 20 0 17 AUTH AC A4 E1 1\n"
        "HETATM '?' O L5 . LIG LA 7 ? 4 0 0 1 20 0 17 AUTH AC A5 E1 1\n#\n"};
    auto reader = mmcif::MmcifReader{input};
    const auto record = reader.next();
    REQUIRE(record.has_value());
    REQUIRE(record->import_metadata.has_value());
    const auto& mapping = *record->import_metadata;
    REQUIRE(mapping.atoms.size() == 5);
    CHECK(mapping.atoms[0].id == "Csite");
    CHECK(mapping.atoms[1].id == "001");
    CHECK(mapping.atoms[2].id == "9223372036854775808");
    CHECK(mapping.atoms[3].id == ".");
    CHECK(mapping.atoms[4].id == "?");
    REQUIRE(mapping.atoms[0].structural_labels.has_value());
    CHECK(mapping.atoms[0].structural_labels->author.atom == "A1");
    CHECK(mapping.atoms[0].structural_labels->author.chain == "AC");
    CHECK(mapping.atoms[0].structural_labels->label.atom == "L1");
    CHECK(mapping.atoms[0].structural_labels->label.chain == "LA");
    CHECK(mapping.atoms[1].structural_labels->alternate_location == "B");

    for (const auto missing : {".", "?"}) {
        std::istringstream missing_input{"data_missing\nloop_\n_atom_site.id\n" +
                                         std::string{missing} + "\n"};
        auto missing_reader = mmcif::MmcifReader{missing_input};
        CHECK_THROWS_AS(missing_reader.next(), std::runtime_error);
    }

    std::istringstream duplicate{
        std::string{header} +
        "HETATM duplicate C C1 . LIG A 1 ? 0 0 0 1 20 0 1 LIG A C1 E1 1\n"
        "HETATM 'duplicate' O O1 . LIG A 1 ? 1 0 0 1 20 0 1 LIG A O1 E1 1\n#\n"};
    auto duplicate_reader = mmcif::MmcifReader{duplicate};
    CHECK_THROWS_AS(duplicate_reader.next(), std::runtime_error);

    std::istringstream long_value{atom_site_header("long") +
                                  "\nHETATM 1 C C1 . LIG A 1 ? 0 0 0 1 20 0 1 LIG A C1 1\n#\n"
                                  "_struct.title " +
                                  std::string(5000, 'x') + "\n"};
    auto long_reader = mmcif::MmcifReader{long_value};
    const auto long_record = long_reader.next();
    REQUIRE(long_record.has_value());
    CHECK(long_record->molecule.atoms().size() == 1);

    for (const auto model_id : {"+1", ".", "?"}) {
        std::istringstream model_input{
            std::string{header} + "HETATM source C L1 . LIG LA 7 ? 0 0 0 1 20 0 17 AUTH AC A1 E1 " +
            model_id + "\n#\n"};
        auto model_reader = mmcif::MmcifReader{model_input};
        const auto model_record = model_reader.next();
        REQUIRE(model_record.has_value());
        REQUIRE(model_record->import_metadata.has_value());
        REQUIRE(model_record->import_metadata->conformers.size() == 1);
        CHECK(model_record->import_metadata->conformers[0].id == model_id);
    }
}

TEST_CASE("mmCIF templates use canonical labels and reject ambiguous selected names",
          "[adapters][mmcif]") {
    const auto make_input = [](const std::string_view atom_rows) {
        return atom_site_header("canonical") + "\n" + std::string{atom_rows} + "\n#\n";
    };

    std::istringstream input{make_input("ATOM 1 N N . ALA A 1 ? 0 0 0 1 20 0 9 AUTHOR A N1 1\n"
                                        "ATOM 2 C CA . ALA A 1 ? 1 0 0 1 20 0 9 AUTHOR A CA1 1")};
    auto reader =
        mmcif::MmcifReader{input, {}, {.bond_strategy = gemmi_adapter::BondStrategy::templates}};
    const auto record = reader.next();
    REQUIRE(record.has_value());
    CHECK(record->molecule.atom(0).name() == "N1");
    CHECK(record->molecule.atom(1).name() == "CA1");
    REQUIRE(record->molecule.bond_count() == 1);
    CHECK(record->molecule.bond(0).order() == chargefw::core::BondOrder::SINGLE);
    REQUIRE(record->import_metadata.has_value());
    CHECK(record->import_metadata->components[0].component_id == "ALA");
    CHECK(record->import_metadata->atoms[0].structural_labels->label.atom == "N");
    CHECK(record->import_metadata->atoms[0].structural_labels->author.atom == "N1");

    std::istringstream terminal{
        make_input("ATOM 1 H H3 . PRO A 1 ? 0 0 0 1 20 0 9 AUTHOR A THIRD 1\n"
                   "ATOM 2 N N . PRO A 1 ? 1 0 0 1 20 0 9 AUTHOR A NITROGEN 1\n"
                   "ATOM 3 H H2 . PRO A 1 ? 2 0 0 1 20 0 9 AUTHOR A SECOND 1\n"
                   "ATOM 4 C CA . PRO A 1 ? 3 0 0 1 20 0 9 AUTHOR A ALPHA 1\n"
                   "ATOM 5 C CD . PRO A 1 ? 4 0 0 1 20 0 9 AUTHOR A DELTA 1")};
    auto terminal_reader =
        mmcif::MmcifReader{terminal, {}, {.bond_strategy = gemmi_adapter::BondStrategy::templates}};
    const auto terminal_record = terminal_reader.next();
    REQUIRE(terminal_record.has_value());
    CHECK(terminal_record->molecule.atom_count() == 5);
    CHECK(terminal_record->molecule.bond_count() == 4);
    for (const auto& bond : terminal_record->molecule.bonds()) {
        CHECK(bond.first_atom_index() == 1);
        CHECK(bond.order() == chargefw::core::BondOrder::SINGLE);
    }
    CHECK(terminal_record->molecule.atom(0).name() == "THIRD");
    CHECK(terminal_record->molecule.atom(1).name() == "NITROGEN");
    CHECK(terminal_record->molecule.atom(1).formal_charge() == 0);
    CHECK(terminal_record->import_metadata->components[0].component_id == "PRO");

    auto explicit_input = make_input("ATOM 1 N N . ALA A 1 ? 0 0 0 1 20 0 9 AUTHOR A N1 1\n"
                                     "ATOM 2 C CA . ALA A 1 ? 1 0 0 1 20 0 9 AUTHOR A CA1 1");
    explicit_input += R"cif(loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
ALA N CA SING
#
)cif";
    std::istringstream explicit_stream{explicit_input};
    auto explicit_reader = mmcif::MmcifReader{
        explicit_stream, {}, {.bond_strategy = gemmi_adapter::BondStrategy::explicit_bonds}};
    const auto explicit_record = explicit_reader.next();
    REQUIRE(explicit_record.has_value());
    REQUIRE(explicit_record->molecule.bond_count() == 1);
    CHECK(explicit_record->molecule.bond(0).order() == chargefw::core::BondOrder::SINGLE);

    std::istringstream fallback{make_input("ATOM 1 N . . . A 1 ? 0 0 0 1 20 0 9 ALA A N 1\n"
                                           "ATOM 2 C . . . A 1 ? 1 0 0 1 20 0 9 ALA A CA 1")};
    auto fallback_reader =
        mmcif::MmcifReader{fallback, {}, {.bond_strategy = gemmi_adapter::BondStrategy::templates}};
    const auto fallback_record = fallback_reader.next();
    REQUIRE(fallback_record.has_value());
    CHECK(fallback_record->molecule.bond_count() == 1);
    CHECK(fallback_record->import_metadata->components[0].component_id == "ALA");

    std::istringstream ambiguous{make_input("HETATM 1 C C1 . LIG A 1 ? 0 0 0 1 20 0 9 LIG A X1 1\n"
                                            "HETATM 2 O C1 . LIG A 1 ? 1 0 0 1 20 0 9 LIG A X2 1")};
    auto ambiguous_reader = mmcif::MmcifReader{ambiguous};
    CHECK_THROWS_AS(ambiguous_reader.next(), std::runtime_error);

    std::istringstream alternates{
        make_input("HETATM 1 C C1 A LIG A 1 ? 0 0 0 1 20 0 9 LIG A X1 1\n"
                   "HETATM 2 C C1 B LIG A 1 ? 1 0 0 1 20 0 9 LIG A X1 1")};
    auto alternate_reader = mmcif::MmcifReader{alternates};
    const auto alternate_record = alternate_reader.next();
    REQUIRE(alternate_record.has_value());
    CHECK(alternate_record->molecule.atom_count() == 1);
}

TEST_CASE("mmCIF input preserves records, models, selection, and bond strategy",
          "[adapters][mmcif]") {
    std::istringstream input{atom_site_header("first") + R"cif(
ATOM 1 C CA . ALA A 1 ? 0.0 0.0 0.0 1.0 20.0 0 1 ALA A CA 1
HETATM 2 C C1 . LIG A 2 ? 1.0 0.0 0.0 1.0 20.0 0 2 LIG A C1 1
HETATM 3 O O . HOH A 3 ? 2.0 0.0 0.0 1.0 20.0 0 3 HOH A O 1
ATOM 4 C CA . ALA A 1 ? 0.1 0.0 0.0 1.0 20.0 0 1 ALA A CA 2
HETATM 5 C C1 . LIG A 2 ? 1.1 0.0 0.0 1.0 20.0 0 2 LIG A C1 2
HETATM 6 O O . HOH A 3 ? 2.1 0.0 0.0 1.0 20.0 0 3 HOH A O 2
#
)cif" + atom_site_header("second") +
                             R"cif(
ATOM 1 O O . HOH B 1 ? 3.0 0.0 0.0 1.0 20.0 0 1 HOH B O 1
#
loop_
_sb_ncbr_partial_atomic_charges.type_id
_sb_ncbr_partial_atomic_charges.atom_id
_sb_ncbr_partial_atomic_charges.charge
1 1 -0.8
#
)cif"};

    auto reader = mmcif::MmcifReader{input, "structure.cif"};
    const auto first = reader.next();
    REQUIRE(first.has_value());
    CHECK(first->identity.source == "structure.cif");
    CHECK(first->identity.record_index == 0);
    CHECK(first->identity.record_id == "first");
    CHECK(first->molecule.atom_count() == 3);
    CHECK(first->molecule.bond_count() == 0);
    CHECK(first->molecule.conformer_count() == 2);
    CHECK(first->molecule.conformer(1)[0].x == 0.1);
    CHECK(first->diagnostics.empty());

    const auto second = reader.next();
    REQUIRE(second.has_value());
    CHECK(second->identity.record_index == 1);
    CHECK(second->identity.record_id == "second");
    CHECK(second->molecule.atom_count() == 1);
    // Charges in the file are not calculation input.
    REQUIRE(second->diagnostics.size() == 1);
    CHECK(second->diagnostics.front().code == "partial_charges_ignored");
    CHECK_FALSE(reader.next().has_value());

    std::istringstream deferred_error_input{atom_site_header("valid") + R"cif(
ATOM 1 C CA . ALA A 1 ? 0.0 0.0 0.0 1.0 20.0 0 1 ALA A CA 1
#
)cif" + atom_site_header("invalid") + R"cif(
ATOM 1 Xx CA . ALA A 1 ? 0.0 0.0 0.0 1.0 20.0 0 1 ALA A CA 1
#
)cif"};
    auto deferred_error_reader = mmcif::MmcifReader{deferred_error_input};
    REQUIRE(deferred_error_reader.next().has_value());
    CHECK_THROWS_AS(deferred_error_reader.next(), std::runtime_error);

    std::istringstream filtered_input{atom_site_header("filtered") + R"cif(
ATOM 1 C CA . ALA A 1 ? 0.0 0.0 0.0 1.0 20.0 0 1 ALA A CA 1
HETATM 2 C C1 . LIG A 2 ? 1.0 0.0 0.0 1.0 20.0 0 2 LIG A C1 1
HETATM 3 O O . HOH A 3 ? 2.0 0.0 0.0 1.0 20.0 0 3 HOH A O 1
#
)cif"};
    auto filtered = mmcif::MmcifReader{
        filtered_input, {}, {.selection = gemmi_adapter::RecordSelection::polymers}};
    const auto filtered_record = filtered.next();
    REQUIRE(filtered_record.has_value());
    CHECK(filtered_record->molecule.atom_count() == 2);

    const auto strategy_input = atom_site_header("connectivity") + R"cif(
HETATM 1 C C1 . LIG A 1 ? 0.0 0.0 0.0 1.0 20.0 0 1 LIG A C1 1
HETATM 2 O O1 . LIG A 1 ? 1.0 0.0 0.0 1.0 20.0 0 1 LIG A O1 1
HETATM 3 C C2 . LIG B 1 ? 2.0 0.0 0.0 1.0 20.0 0 1 LIG B C2 1
HETATM 4 O O2 . LIG B 1 ? 3.0 0.0 0.0 1.0 20.0 0 1 LIG B O2 1
#
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
LIG C1 O1 DOUB
LIG C2 O2 SING
#
loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_comp_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_comp_id
_struct_conn.ptnr2_label_seq_id
_struct_conn.ptnr2_label_atom_id
_struct_conn.pdbx_value_order
link1 covale A LIG 1 O1 B LIG 1 C2 doub
link2 covale A LIG 1 C1 B LIG 1 O2 ?
#
)cif";
    const auto read_strategy = [&](const gemmi_adapter::BondStrategy strategy) {
        return read_first(strategy_input, {.bond_strategy = strategy});
    };

    CHECK(read_strategy(gemmi_adapter::BondStrategy::none).bond_count() == 0);
    const auto explicit_molecule = read_strategy(gemmi_adapter::BondStrategy::explicit_bonds);
    REQUIRE(explicit_molecule.bond_count() == 4);
    CHECK(explicit_molecule.bond(0).order() == chargefw::core::BondOrder::DOUBLE);
    CHECK(explicit_molecule.bond(1).order() == chargefw::core::BondOrder::SINGLE);
    // Connections between residues take their pdbx_value_order, single when it is missing.
    CHECK(explicit_molecule.bond(2).order() == chargefw::core::BondOrder::DOUBLE);
    CHECK(explicit_molecule.bond(3).order() == chargefw::core::BondOrder::SINGLE);
    CHECK(read_strategy(gemmi_adapter::BondStrategy::hybrid).bond_count() == 4);

    // A connection to an atom of a symmetry-related copy is not a bond of this molecule.
    const auto read_connection = [](const std::string_view second_symmetry) {
        std::istringstream connection_input{atom_site_header("symmetry") + R"cif(
HETATM 1 S S1 . LIG A 1 ? 0.0 0.0 0.0 1.0 20.0 0 1 LIG A S1 1
HETATM 2 S S2 . LIG B 1 ? 15.0 0.0 0.0 1.0 20.0 0 1 LIG B S2 1
#
loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_comp_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr1_symmetry
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_comp_id
_struct_conn.ptnr2_label_seq_id
_struct_conn.ptnr2_label_atom_id
_struct_conn.ptnr2_symmetry
disulf1 disulf A LIG 1 S1 1_555 B LIG 1 S2 )cif" +
                                            std::string{second_symmetry} + "\n#\n"};
        auto record =
            mmcif::MmcifReader{connection_input,
                               {},
                               {.bond_strategy = gemmi_adapter::BondStrategy::explicit_bonds}}
                .next();
        REQUIRE(record.has_value());
        return std::move(*record);
    };
    const auto same_copy = read_connection("1_555");
    CHECK(same_copy.molecule.bond_count() == 1);
    CHECK(same_copy.diagnostics.empty());
    const auto other_copy = read_connection("2_555");
    CHECK(other_copy.molecule.bond_count() == 0);
    REQUIRE(other_copy.diagnostics.size() == 1);
    CHECK(other_copy.diagnostics.front().code == "symmetry_connections_ignored");

    const auto duplicate_input = [](const std::string_view component_bond) {
        return atom_site_header("duplicate") + R"cif(
ATOM 1 N N . ALA A 1 ? 0.0 0.0 0.0 1.0 20.0 0 1 ALA A N 1
ATOM 2 C CA . ALA A 1 ? 1.0 0.0 0.0 1.0 20.0 0 1 ALA A CA 1
#
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
)cif" + std::string{component_bond} +
               "\n#\n";
    };
    const auto read_duplicate_bond = [&](const gemmi_adapter::BondStrategy strategy,
                                         const std::string_view component_bond) {
        const auto molecule =
            read_first(duplicate_input(component_bond), {.bond_strategy = strategy});
        REQUIRE(molecule.bond_count() == 1);
        return molecule.bond(0);
    };

    CHECK(read_duplicate_bond(gemmi_adapter::BondStrategy::templates, "ALA N CA DOUB").order() ==
          chargefw::core::BondOrder::SINGLE);
    for (const auto strategy :
         {gemmi_adapter::BondStrategy::explicit_bonds, gemmi_adapter::BondStrategy::hybrid}) {
        const auto unquoted = read_duplicate_bond(strategy, "ALA N CA DOUB");
        const auto quoted = read_duplicate_bond(strategy, "'ALA' 'N' 'CA' 'DOUB'");
        CHECK(quoted.first_atom_index() == unquoted.first_atom_index());
        CHECK(quoted.second_atom_index() == unquoted.second_atom_index());
        CHECK(quoted.order() == unquoted.order());
        CHECK(quoted.order() == chargefw::core::BondOrder::DOUBLE);
    }

    const auto explicit_bond_count = [&](const std::string_view value_order) {
        return read_first(duplicate_input("ALA N CA " + std::string{value_order}),
                          {.bond_strategy = gemmi_adapter::BondStrategy::explicit_bonds})
            .bond_count();
    };
    CHECK(
        read_duplicate_bond(gemmi_adapter::BondStrategy::explicit_bonds, "ALA N CA aRoM").order() ==
        chargefw::core::BondOrder::SINGLE);
    for (const auto value_order : {"DELO", "METAL", "1.5", ".", "?", "invalid"}) {
        CHECK(explicit_bond_count(value_order) == 0);
    }
}

TEST_CASE("mmCIF input restores source order across residues and models", "[adapters][mmcif]") {
    std::istringstream input{atom_site_header("order") + R"cif(
HETATM Csite C C1 . LIG A 1 ? 0 0 0 1 20 0 1 LIG A C1 1
HETATM Osite O O1 . LIG A 2 ? 1 0 0 1 20 0 2 LIG A O1 1
HETATM Nsite N N1 . LIG A 1 ? 2 0 0 1 20 0 1 LIG A N1 1
HETATM Csite2 C C1 . LIG A 1 ? 10 0 0 1 20 0 1 LIG A C1 2
HETATM Osite2 O O1 . LIG A 2 ? 11 0 0 1 20 0 2 LIG A O1 2
HETATM Nsite2 N N1 . LIG A 1 ? 12 0 0 1 20 0 1 LIG A N1 2
#
loop_
_chem_comp_bond.comp_id
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
LIG C1 N1 SING
#
)cif"};
    auto reader = mmcif::MmcifReader{
        input, {}, {.bond_strategy = gemmi_adapter::BondStrategy::explicit_bonds}};
    const auto record = reader.next();
    REQUIRE(record.has_value());
    REQUIRE(record->molecule.atom_count() == 3);
    CHECK(record->molecule.atom(0).name() == "C1");
    CHECK(record->molecule.atom(1).name() == "O1");
    CHECK(record->molecule.atom(2).name() == "N1");
    CHECK(record->molecule.conformer(0)[1].x == 1.0);
    CHECK(record->molecule.conformer(1)[0].x == 10.0);
    CHECK(record->molecule.conformer(1)[1].x == 11.0);
    CHECK(record->molecule.conformer(1)[2].x == 12.0);
    REQUIRE(record->molecule.bond_count() == 1);
    CHECK(record->molecule.bond(0).first_atom_index() == 0);
    CHECK(record->molecule.bond(0).second_atom_index() == 2);
    REQUIRE(record->import_metadata.has_value());
    CHECK(record->import_metadata->source_connectivity ==
          chargefw::adapters::SourceConnectivity::present);
    CHECK(record->import_metadata->atoms[0].id == "Csite");
    CHECK(record->import_metadata->atoms[1].id == "Osite");
    CHECK(record->import_metadata->atoms[2].id == "Nsite");
    REQUIRE(record->import_metadata->components.size() == 2);
    CHECK(record->import_metadata->components[0] ==
          chargefw::adapters::SourceComponentInstance{"LIG", {0, 2}});
    CHECK(record->import_metadata->components[1] ==
          chargefw::adapters::SourceComponentInstance{"LIG", {1}});
}

TEST_CASE("mmCIF input rejects incompatible conformer atom sequences", "[adapters][mmcif]") {
    const auto incompatible_input = atom_site_header("incompatible") + R"cif(
HETATM 1 C C1 . LIG A 1 ? 0.0 0.0 0.0 1.0 20.0 0 1 LIG A C1 1
HETATM 2 O O1 . LIG A 1 ? 1.0 0.0 0.0 1.0 20.0 0 1 LIG A O1 1
HETATM 3 C C1 . LIG A 1 ? 0.1 0.0 0.0 1.0 20.0 0 1 LIG A C1 2
HETATM 4 O O1 . LIG A 1 ? 1.1 0.0 0.0 1.0 20.0 0 1 LIG B O1 2
#
)cif";
    std::istringstream input{incompatible_input};
    auto reader = mmcif::MmcifReader{input};
    CHECK_THROWS_MATCHES(reader.next(), std::exception,
                         snitch::matchers::with_what_contains{"same selected component instances"});
}
