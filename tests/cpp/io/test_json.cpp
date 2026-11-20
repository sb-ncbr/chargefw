#include "support/test_exceptions.h"

#include <chargefw/core/bond.h>
#include <chargefw/io/json_input.h>
#include <chargefw/io/parse_error.h>
#include <snitch/snitch.hpp>

#include <array>
#include <exception>
#include <sstream>
#include <string>
#include <type_traits>

namespace io = chargefw::io;

static_assert(!std::is_copy_constructible_v<io::JsonReader> &&
              !std::is_copy_assignable_v<io::JsonReader>);
static_assert(std::is_move_constructible_v<io::JsonReader>);

TEST_CASE("JSON input preserves identity, graph, and conformer mapping", "[io][json]") {
    {
        std::istringstream input{R"json(
{
  "schema_version": "1.0",
  "molecules": [
    {
      "id": "water-1",
      "name": "water",
      "atoms": [
        { "atomic_number": 8, "formal_charge": 0 },
        { "atomic_number": 1, "formal_charge": 0 },
        { "atomic_number": 1, "formal_charge": 0 }
      ],
      "bonds": [
        { "atoms": [0, 1], "order": 1 },
        { "atoms": [0, 2], "order": 1 }
      ],
      "conformers": [
        {
          "id": "model-1",
          "coordinates": [
            [0.0, 0.0, 0.0],
            [0.9572, 0.0, 0.0],
            [-0.239987, 0.927297, 0.0]
          ]
        }
      ]
    }
  ]
}
)json"};
        auto reader = io::JsonReader{input, "water.json"};
        const auto result = reader.next();
        REQUIRE(result.has_value());
        const auto& record = *result;
        CHECK(record.identity.source_name == "water.json");
        CHECK(record.identity.record_index == 0);
        CHECK(record.identity.record_id == "water-1");
        CHECK(record.molecule.name() == "water");
        CHECK(record.molecule.atom_count() == 3);
        CHECK(record.molecule.atom(0).atomic_number() == 8);
        CHECK(record.molecule.atom(0).formal_charge() == 0);
        CHECK(record.molecule.atom(0).name().empty());
        CHECK(record.molecule.bond_count() == 2);
        CHECK(record.molecule.bond(0).order() == chargefw::core::BondOrder::SINGLE);
        CHECK(record.molecule.conformer_count() == 1);
        CHECK(record.molecule.conformer(0).name() == "model-1");
        CHECK(record.molecule.conformer(0)[1].x == 0.9572);
        CHECK_FALSE(reader.next().has_value());
    }

    {
        std::istringstream input{R"json(
{
  "schema_version": "1.0",
  "molecules": [
    {
      "id": "valid",
      "atoms": [{ "atomic_number": 1, "formal_charge": 0 }]
    },
    {
      "id": "invalid",
      "atoms": [{ "atomic_number": 6, "formal_charge": 0 }],
      "bonds": [{ "atoms": [0, 1], "order": 1 }]
    }
  ]
}
)json"};
        auto reader = io::JsonReader{input, "input.json"};
        CHECK(reader.next().has_value());
        // The error gives the record; JSON has no lines.
        const auto error =
            chargefw::test::thrown<io::ParseError>([&] { static_cast<void>(reader.next()); });
        REQUIRE(error.has_value());
        CHECK(error->location().source_name == "input.json");
        CHECK(error->location().record_index == 1);
        CHECK(error->location().record_id == "invalid");
        CHECK_FALSE(error->location().line.has_value());
    }

    {
        std::istringstream input{R"json(
{
  "schema_version": "2.0",
  "molecules": []
}
)json"};
        CHECK_THROWS_AS(io::JsonReader{input}, std::runtime_error);
    }
}

TEST_CASE("JSON input tells syntax and version apart", "[io][json]") {
    const auto kind = [](const std::string& text) {
        std::istringstream input{text};
        const auto error = chargefw::test::thrown<io::ParseError>(
            [&] { static_cast<void>(io::JsonReader{input, "input.json"}); });
        REQUIRE(error.has_value());
        CHECK(error->location().source_name == "input.json");
        return error->kind();
    };
    // JSON is UTF-8, so a Latin-1 byte breaks the format as a missing comma does.
    CHECK(kind("{\"schema_version\": \"1.0\", \"molecules\": [{\"id\": \"caf\xe9\"}]}") ==
          io::ParseErrorKind::malformed);
    CHECK(kind("{\"schema_version\": \"1.0\" \"molecules\": [{\"id\": \"caf\xc3\xa9\"}]}") ==
          io::ParseErrorKind::malformed);
    // Which members another version allows is not known, so the version is what is refused.
    CHECK(kind(R"json({"schema_version": "2.0", "molecules": [], "units": "angstrom"})json") ==
          io::ParseErrorKind::unsupported);
    CHECK(kind(R"json({"schema_version": "1.0", "molecules": [], "units": "angstrom"})json") ==
          io::ParseErrorKind::malformed);
}

TEST_CASE("JSON input passes a failed read on instead of reading or refusing the input",
          "[io][json]") {
    for (const auto state : {std::ios::failbit, std::ios::badbit}) {
        std::istringstream input{R"json({"schema_version": "1.0", "molecules": []})json"};
        input.setstate(state);
        auto read = false;
        auto refused = false;
        auto failed = false;
        try {
            static_cast<void>(io::JsonReader{input, "input.json"});
            read = true;
        } catch (const io::ParseError&) {
            refused = true;
        } catch (const std::runtime_error&) {
            failed = true;
        }
        CHECK(failed);
        CHECK_FALSE(read);
        CHECK_FALSE(refused);
    }
}

TEST_CASE("JSON input rejects integers outside the native integer range", "[io][json]") {
    const auto input = R"json(
{
  "schema_version": "1.0",
  "molecules": [
    {"atoms": [{"atomic_number": 4294967297, "formal_charge": 0}]}
  ]
}
)json";
    std::istringstream stream{input};
    auto reader = io::JsonReader{stream};
    CHECK_THROWS_AS(reader.next(), std::runtime_error);
}

TEST_CASE("JSON input rejects unknown members", "[io][json]") {
    // One misspelled or unsupported member per object kind.
    const std::array documents{
        R"json({"schema_version": "1.0", "molecules": [], "molecule": []})json",
        R"json({"schema_version": "1.0", "molecules": [
            {"atoms": [{"atomic_number": 1, "formal_charge": 0}], "bond": []}]})json",
        R"json({"schema_version": "1.0", "molecules": [
            {"atoms": [{"atomic_number": 1, "formal_charge": 0, "charge": 0.5}]}]})json",
        R"json({"schema_version": "1.0", "molecules": [
            {"atoms": [{"atomic_number": 1, "formal_charge": 0},
                       {"atomic_number": 1, "formal_charge": 0}],
             "bonds": [{"atoms": [0, 1], "order": 1, "aromatic": false}]}]})json",
        R"json({"schema_version": "1.0", "molecules": [
            {"atoms": [{"atomic_number": 1, "formal_charge": 0}],
             "conformers": [{"coordinates": [[0, 0, 0]], "units": "nm"}]}]})json",
    };
    for (const auto* document : documents) {
        CAPTURE(document);
        std::istringstream input{document};
        CHECK_THROWS_AS(static_cast<void>(io::JsonReader{input}.next()), std::runtime_error);
    }

    // Also in a conformer that reading only the first one skips.
    std::istringstream skipped{R"json({"schema_version": "1.0", "molecules": [
        {"atoms": [{"atomic_number": 1, "formal_charge": 0}],
         "conformers": [{"coordinates": [[0, 0, 0]]},
                        {"coordinates": [[0, 0, 0]], "units": "nm"}]}]})json"};
    CHECK_THROWS_AS(
        static_cast<void>(
            io::JsonReader{skipped, {}, chargefw::io::ConformerSelection::first}.next()),
        std::runtime_error);
}

TEST_CASE("JSON input owns coordinate-free source mapping", "[io][json]") {
    auto record = [] {
        std::istringstream input{R"json(
{
  "schema_version": "1.0",
  "molecules": [{
    "atoms": [
      {"atomic_number": 6, "formal_charge": 0},
      {"atomic_number": 8, "formal_charge": 0}
    ]
  }]
}
)json"};
        auto reader = io::JsonReader{input, "topology.json"};
        return *reader.next();
    }();

    REQUIRE(record.import_metadata.has_value());
    const auto& metadata = *record.import_metadata;
    CHECK(metadata.format == chargefw::io::MolecularSourceFormat::molecule_json);
    REQUIRE(metadata.atoms.size() == 2);
    CHECK(metadata.atoms[0].position == 0);
    CHECK_FALSE(metadata.atoms[0].id.has_value());
    CHECK(metadata.atoms[1].position == 1);
    CHECK(metadata.conformers.empty());
    CHECK(metadata.policy.conformers == "all");
    CHECK(metadata.source_connectivity == chargefw::io::SourceConnectivity::absent);
}
