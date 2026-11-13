#include <chargefw/adapters/native/json_input.h>
#include <chargefw/core/bond.h>
#include <snitch/snitch.hpp>

#include <array>
#include <exception>
#include <sstream>
#include <string>
#include <type_traits>

namespace adapters = chargefw::adapters;

static_assert(!std::is_copy_constructible_v<adapters::JsonReader> &&
              !std::is_copy_assignable_v<adapters::JsonReader>);
static_assert(std::is_move_constructible_v<adapters::JsonReader>);

TEST_CASE("JSON input preserves identity, graph, and conformer mapping", "[adapters][json]") {
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
        auto reader = adapters::JsonReader{input, "water.json"};
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
      "id": "invalid",
      "atoms": [{ "atomic_number": 6, "formal_charge": 0 }],
      "bonds": [{ "atoms": [0, 1], "order": 1 }]
    },
    {
      "id": "valid",
      "atoms": [{ "atomic_number": 1, "formal_charge": 0 }]
    }
  ]
}
)json"};
        auto reader = adapters::JsonReader{input};
        CHECK_THROWS_AS(reader.next(), std::exception);
    }

    {
        std::istringstream input{R"json(
{
  "schema_version": "2.0",
  "molecules": []
}
)json"};
        CHECK_THROWS_AS(adapters::JsonReader{input}, std::runtime_error);
    }
}

TEST_CASE("JSON input rejects integers outside the native integer range", "[adapters][json]") {
    const auto input = R"json(
{
  "schema_version": "1.0",
  "molecules": [
    {"atoms": [{"atomic_number": 4294967297, "formal_charge": 0}]}
  ]
}
)json";
    std::istringstream stream{input};
    auto reader = adapters::JsonReader{stream};
    CHECK_THROWS_AS(reader.next(), std::runtime_error);
}

TEST_CASE("JSON input rejects unknown members", "[adapters][json]") {
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
        CHECK_THROWS_AS(static_cast<void>(adapters::JsonReader{input}.next()), std::runtime_error);
    }

    // Also in a conformer that reading only the first one skips.
    std::istringstream skipped{R"json({"schema_version": "1.0", "molecules": [
        {"atoms": [{"atomic_number": 1, "formal_charge": 0}],
         "conformers": [{"coordinates": [[0, 0, 0]]},
                        {"coordinates": [[0, 0, 0]], "units": "nm"}]}]})json"};
    CHECK_THROWS_AS(static_cast<void>(adapters::JsonReader{
                        skipped, {}, chargefw::adapters::ConformerSelection::first}
                                          .next()),
                    std::runtime_error);
}

TEST_CASE("JSON input owns coordinate-free source mapping", "[adapters][json]") {
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
        auto reader = adapters::JsonReader{input, "topology.json"};
        return *reader.next();
    }();

    REQUIRE(record.import_metadata.has_value());
    const auto& metadata = *record.import_metadata;
    CHECK(metadata.format == chargefw::adapters::MolecularSourceFormat::molecule_json);
    REQUIRE(metadata.atoms.size() == 2);
    CHECK(metadata.atoms[0].position == 0);
    CHECK_FALSE(metadata.atoms[0].id.has_value());
    CHECK(metadata.atoms[1].position == 1);
    CHECK(metadata.conformers.empty());
    CHECK(metadata.policy.conformers == "all");
    CHECK(metadata.source_connectivity == chargefw::adapters::SourceConnectivity::absent);
}
