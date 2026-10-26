#include <chargefw/methods/method_registry.h>
#include <chargefw/parameters/io/parameter_set_io.h>

#include <filesystem>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_set>

#include <snitch/snitch.hpp>

namespace parameters = chargefw::parameters;

namespace {

#ifndef CHARGEFW_TEST_PARAMETER_DIR
#error "CHARGEFW_TEST_PARAMETER_DIR must be defined"
#endif

auto peoe_file() -> std::filesystem::path {
    return std::filesystem::path{CHARGEFW_TEST_PARAMETER_DIR} / "PEOE_original.json";
}

auto bad_length_json() -> std::string {
    return R"json(
{
  "metadata": {
    "id": "bad",
    "name": "Bad",
    "method": "peoe"
  },
  "atom": {
    "names": ["A", "B", "C"],
    "data": [
      {
        "key": {"element": "H", "classifier": "hbo", "type": "1"},
        "values": [7.17, 6.24]
      }
    ]
  }
}
)json";
}

auto missing_id_json() -> std::string {
    return R"json(
{
  "metadata": {
    "name": "No ID",
    "method": "peoe"
  }
}
)json";
}

auto named_bond_key_json() -> std::string {
    return R"json(
{
  "metadata": {
    "id": "named-bond",
    "name": "Named bond",
    "method": "sqeqp"
  },
  "bond": {
    "names": ["kappa"],
    "data": [
      {
        "key": {
          "atoms": [
            {"element": "C", "classifier": "bonded", "type": "CNO"},
            {"element": "C", "classifier": "bonded", "type": "CCN"}
          ],
          "bond": {"classifier": "bo", "type": "1"}
        },
        "values": [1.0]
      }
    ]
  }
}
)json";
}

auto unknown_element_json() -> std::string {
    return R"json(
{
  "metadata": {
    "id": "bad",
    "name": "Bad",
    "method": "peoe"
  },
  "atom": {
    "names": ["A"],
    "data": [
      {
        "key": {"element": "Xx", "classifier": "hbo", "type": "1"},
        "values": [1.0]
      }
    ]
  }
}
)json";
}

auto unknown_classification_json() -> std::string {
    return R"json(
{
  "metadata": {
    "id": "bad",
    "name": "Bad",
    "method": "peoe"
  },
  "atom": {
    "names": ["A"],
    "data": [
      {
        "key": {"element": "H", "classifier": "unknown", "type": "1"},
        "values": [1.0]
      }
    ]
  }
}
)json";
}

} // namespace

TEST_CASE("parameter set loads from a PEOE JSON file", "[parameters][io]") {
    const auto parameter_set = parameters::load_parameter_set_json_file(peoe_file());

    CHECK(parameter_set.id() == std::string_view{"PEOE_original"});
    CHECK(parameter_set.method_id() == std::string_view{"peoe"});

    CHECK(parameter_set.common().size() == 1);
    CHECK(parameter_set.common().contains("dampH"));
    CHECK(parameter_set.common().parameter("dampH") == 20.02);

    REQUIRE_FALSE(parameter_set.atom().empty());

    const auto& hydrogen = parameter_set.atom()[0];

    CHECK(hydrogen.key.atomic_number == 1);
    CHECK(hydrogen.key.classification ==
          parameters::AtomParameterClassificationKind::HIGHEST_BOND_ORDER);
    CHECK(hydrogen.key.type == "1");

    CHECK(parameter_set.atom().contains(0, "A"));
    CHECK(parameter_set.atom().contains(0, "B"));
    CHECK(parameter_set.atom().contains(0, "C"));

    CHECK(parameter_set.atom().parameter(0, "A") == 7.17);
    CHECK(parameter_set.atom().parameter(0, "B") == 6.24);
    CHECK(parameter_set.atom().parameter(0, "C") == -0.56);
}

TEST_CASE("parameter sets load grouped by method in selection order", "[parameters][io]") {
    const auto parameter_sets =
        parameters::load_parameter_sets_json_directory(CHARGEFW_TEST_PARAMETER_DIR);
    for (std::size_t index = 1; index < parameter_sets.size(); ++index) {
        const auto& previous = parameter_sets[index - 1];
        const auto& current = parameter_sets[index];
        CAPTURE(previous.id(), current.id());
        CHECK(previous.method_id() <= current.method_id());
        if (previous.method_id() == current.method_id()) {
            CHECK(previous.priority() >= current.priority());
            if (previous.priority() == current.priority()) {
                CHECK(previous.id() < current.id());
            }
        }
    }
}

TEST_CASE("bundled parameter sets identify registered methods and satisfy their requirements",
          "[parameters][io]") {
    const auto parameter_sets =
        parameters::load_parameter_sets_json_directory(CHARGEFW_TEST_PARAMETER_DIR);
    REQUIRE_FALSE(parameter_sets.empty());
    const auto& registry = chargefw::methods::method_registry();
    auto ids = std::unordered_set<std::string>{};

    for (const auto& parameter_set : parameter_sets) {
        CAPTURE(parameter_set.id());
        CHECK_FALSE(parameter_set.id().empty());
        CHECK(ids.insert(std::string{parameter_set.id()}).second);

        const auto* method = registry.find(parameter_set.method_id());
        REQUIRE(method != nullptr);
        const auto requirements = method->requirements();

        for (const auto name : requirements.common_parameters) {
            CHECK(parameter_set.common().contains(name));
        }

        if (!requirements.atom_parameters.empty()) {
            REQUIRE_FALSE(parameter_set.atom().empty());
        }
        for (std::size_t entry_index = 0; entry_index < parameter_set.atom().size();
             ++entry_index) {
            for (const auto name : requirements.atom_parameters) {
                CHECK(parameter_set.atom().contains(entry_index, name));
            }
        }

        if (!requirements.bond_parameters.empty()) {
            REQUIRE_FALSE(parameter_set.bond().empty());
        }
        for (std::size_t entry_index = 0; entry_index < parameter_set.bond().size();
             ++entry_index) {
            for (const auto name : requirements.bond_parameters) {
                CHECK(parameter_set.bond().contains(entry_index, name));
            }
        }
    }
}

TEST_CASE("parameter set rejects malformed JSON", "[parameters][io]") {
    const auto load_missing_id = [] {
        std::istringstream input{missing_id_json()};
        static_cast<void>(parameters::load_parameter_set_json(input));
    };
    CHECK_THROWS_AS(load_missing_id(), std::invalid_argument);

    const auto load_bad_length = [] {
        std::istringstream input{bad_length_json()};
        static_cast<void>(parameters::load_parameter_set_json(input));
    };
    CHECK_THROWS_AS(load_bad_length(), std::invalid_argument);

    const auto load_unknown_element = [] {
        std::istringstream input{unknown_element_json()};
        static_cast<void>(parameters::load_parameter_set_json(input));
    };
    CHECK_THROWS_AS(load_unknown_element(), std::invalid_argument);

    const auto load_unknown_classification = [] {
        std::istringstream input{unknown_classification_json()};
        static_cast<void>(parameters::load_parameter_set_json(input));
    };
    CHECK_THROWS_MATCHES(load_unknown_classification(), std::invalid_argument,
                         snitch::matchers::with_what_contains{"atom.data[0].key.classifier"});

    const auto load_section = [](const std::string_view section) {
        std::istringstream input{
            R"({"metadata": {"id": "bad", "name": "Bad", "method": "peoe"}, )" +
            std::string{section} + "}"};
        static_cast<void>(parameters::load_parameter_set_json(input));
    };
    CHECK_THROWS_MATCHES(
        load_section(R"("atom": {"names": ["A"], "data": [{"key": )"
                     R"({"element": "H", "classifier": "plain", "type": "1"}, "values": [1.0]}]})"),
        std::invalid_argument, snitch::matchers::with_what_contains{"atom.data[0].key.type"});
    CHECK_THROWS_MATCHES(
        load_section(R"("bond": {"names": ["B"], "data": [{"key": {"atoms": [)"
                     R"({"element": "H", "classifier": "plain", "type": "*"}, )"
                     R"({"element": "H", "classifier": "plain", "type": "*"}], )"
                     R"("bond": {"classifier": "plain", "type": "1"}}, "values": [1.0]}]})"),
        std::invalid_argument, snitch::matchers::with_what_contains{"bond.data[0].key.bond.type"});
}

TEST_CASE("parameter set rejects priorities outside uint16 range", "[parameters][io]") {
    const auto load = [](const std::string_view priority) {
        std::istringstream input{"{\"metadata\":{\"id\":\"test\",\"name\":\"Test\","
                                 "\"method\":\"peoe\",\"priority\":" +
                                 std::string{priority} + "}}"};
        static_cast<void>(parameters::load_parameter_set_json(input));
    };

    CHECK_THROWS_AS(load("-1"), std::invalid_argument);
    CHECK_THROWS_AS(load("65536"), std::invalid_argument);
    CHECK_THROWS_AS(load("4294967297"), std::invalid_argument);
    CHECK_NOTHROW(load("0"));
    CHECK_NOTHROW(load("65535"));
}

TEST_CASE("parameter set loads named unordered bond keys", "[parameters][io]") {
    std::istringstream input{named_bond_key_json()};
    const auto parameter_set = parameters::load_parameter_set_json(input);

    REQUIRE(parameter_set.bond().size() == 1);
    CHECK(parameter_set.bond().parameter(0, "kappa") == 1.0);
}
