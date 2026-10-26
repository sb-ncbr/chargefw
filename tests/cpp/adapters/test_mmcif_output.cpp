#include "support/test_charge_results.h"

#include <chargefw/adapters/charge_result.h>
#include <chargefw/adapters/gemmi/mmcif_input.h>
#include <chargefw/adapters/gemmi/mmcif_output.h>
#include <chargefw/charges/atomic_charges.h>
#include <chargefw/charges/charge_collection.h>

#include <gemmi/cif.hpp>
#include <gemmi/read_cif.hpp>
#include <snitch/snitch.hpp>

#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace {

namespace adapters = chargefw::adapters;
namespace charges = chargefw::charges;
namespace core = chargefw::core;
namespace mmcif_input = chargefw::adapters::gemmi::mmcif_input;
namespace mmcif_output = chargefw::adapters::gemmi::mmcif_output;

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

} // namespace

TEST_CASE("result-owned mmCIF preserves hierarchy and joins conformer charges",
          "[adapters][mmcif]") {
    const auto result = structural_result();
    auto output = std::ostringstream{};
    mmcif_output::MmcifWriter{output}.write(result, "test");
    auto document = ::gemmi::cif::read_string(output.str());
    REQUIRE(document.blocks.size() == 1);
    auto& block = document.blocks.front();
    CHECK_FALSE(block.has_mmcif_category("_chem_comp."));
    CHECK_FALSE(block.has_mmcif_category("_chem_comp_bond."));

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
    REQUIRE(metadata.length() == 2);
    CHECK(::gemmi::cif::as_string(metadata[0][0]) == "ChargeFW");
    CHECK(::gemmi::cif::as_string(metadata[0][1]) == "test");
    CHECK(::gemmi::cif::as_string(values[0][0]) == "1");
    CHECK(::gemmi::cif::as_string(values[2][0]) == "2");
    CHECK(::gemmi::cif::as_string(values[2][1]) == "3");
    CHECK(::gemmi::cif::as_string(values[0][2]) == "0.1235");

    auto round_trip_stream = std::istringstream{output.str()};
    auto round_trip_reader = mmcif_input::MmcifReader{round_trip_stream};
    const auto round_trip = round_trip_reader.next();
    REQUIRE(round_trip.has_value());
    CHECK(round_trip->molecule.atom_count() == 2);
    CHECK(round_trip->molecule.conformer_count() == 2);
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
