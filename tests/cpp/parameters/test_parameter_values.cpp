#include "support/test_parameters.h"

#include <chargefw/parameters/models/atom_parameters.h>
#include <chargefw/parameters/models/bond_parameters.h>
#include <chargefw/parameters/models/common_parameters.h>
#include <chargefw/parameters/models/parameter_set.h>
#include <chargefw/parameters/models/parameter_set_metadata.h>
#include <chargefw/parameters/models/parameter_values.h>

#include <array>
#include <limits>
#include <snitch/snitch.hpp>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace parameters = chargefw::parameters;
using parameters::ParameterScope;

namespace {

struct Variant {
    std::string method_id = "test";
    std::vector<parameters::NamedParameter> common = {{.name = "kappa", .value = 0.5}};
    std::vector<parameters::AtomParameterEntry> atoms = {
        {.key = chargefw::test::plain_atom_key(1),
         .parameters = {{.name = "A", .value = 1.0}, {.name = "B", .value = 2.0}}},
        {.key = chargefw::test::plain_atom_key(8),
         .parameters = {{.name = "A", .value = 3.0}, {.name = "B", .value = 4.0}}}};
    std::vector<parameters::BondParameterEntry> bonds = {
        {.key = chargefw::test::single_bond_key(1, 8),
         .parameters = {{.name = "kappa", .value = 5.0}}}};
};

auto make_set(Variant variant, std::string id = "test-values") -> parameters::ParameterSet {
    return parameters::ParameterSet{
        {.id = std::move(id), .method_id = std::move(variant.method_id)},
        parameters::CommonParameters{std::move(variant.common)},
        parameters::AtomParameters{std::move(variant.atoms)},
        parameters::BondParameters{std::move(variant.bonds)}};
}

} // namespace

TEST_CASE("parameter layout lists common, atom, and bond values in order",
          "[parameters][parameter-values]") {
    const auto set = make_set({});

    const auto expected = std::vector<parameters::ParameterSlot>{
        {.scope = ParameterScope::common, .entry_index = 0, .name = "kappa"},
        {.scope = ParameterScope::atom, .entry_index = 0, .name = "A"},
        {.scope = ParameterScope::atom, .entry_index = 0, .name = "B"},
        {.scope = ParameterScope::atom, .entry_index = 1, .name = "A"},
        {.scope = ParameterScope::atom, .entry_index = 1, .name = "B"},
        {.scope = ParameterScope::bond, .entry_index = 0, .name = "kappa"}};
    CHECK(parameters::parameter_layout(set) == expected);
    CHECK(parameters::parameter_values(set) == std::vector{0.5, 1.0, 2.0, 3.0, 4.0, 5.0});
}

TEST_CASE("new parameter values keep the keys and names of the set",
          "[parameters][parameter-values]") {
    const auto set = make_set({});
    const auto values = std::array{1.5, -1.0, 2.5, -3.0, 4.5, 6.0};

    const auto candidate =
        parameters::with_parameter_values(set, values, {.id = "candidate", .name = "Candidate"});

    CHECK(candidate.id() == "candidate");
    CHECK(candidate.name() == "Candidate");
    CHECK(candidate.method_id() == "test");
    CHECK(parameters::parameter_values(candidate) == std::vector(values.begin(), values.end()));
    CHECK(parameters::parameter_layout(candidate) == parameters::parameter_layout(set));
    CHECK(parameters::compatible(set, candidate));
}

TEST_CASE("parameter sets with other keys, names, order, or method are incompatible",
          "[parameters][parameter-values]") {
    const auto reference = make_set({});
    const auto base = Variant{};

    auto other_method = base;
    other_method.method_id = "other";

    auto swapped_entries = base;
    std::swap(swapped_entries.atoms[0], swapped_entries.atoms[1]);

    auto reordered_names = base;
    std::swap(reordered_names.atoms[0].parameters[0], reordered_names.atoms[0].parameters[1]);

    auto missing_entry = base;
    missing_entry.bonds.clear();

    auto other_common = base;
    other_common.common = {{.name = "sigma", .value = 0.5}};

    // Same entry order, but one key matches other atoms or bonds.
    auto other_atom_classifier = base;
    other_atom_classifier.atoms[1].key = chargefw::test::hbo_atom_key(8, "1");

    auto other_bond_type = base;
    other_bond_type.bonds[0].key.bond.type = "2";

    for (const auto& [label, variant] : std::array{
             std::pair{"other method", other_method}, std::pair{"swapped entries", swapped_entries},
             std::pair{"reordered names", reordered_names},
             std::pair{"missing entry", missing_entry},
             std::pair{"other common parameter", other_common},
             std::pair{"other atom classifier", other_atom_classifier},
             std::pair{"other bond type", other_bond_type}}) {
        CAPTURE(label);
        CHECK_FALSE(parameters::compatible(reference, make_set(variant)));
    }
}

TEST_CASE("new parameter values are rejected for a wrong count, non-finite value, or method",
          "[parameters][parameter-values]") {
    const auto set = make_set({});

    CHECK_THROWS_AS(parameters::with_parameter_values(set, std::array{1.0}, {.id = "short"}),
                    std::invalid_argument);

    auto values = parameters::parameter_values(set);
    values[3] = std::numeric_limits<double>::quiet_NaN();
    CHECK_THROWS_AS(parameters::with_parameter_values(set, values, {.id = "nan"}),
                    std::invalid_argument);

    CHECK_THROWS_AS(parameters::with_parameter_values(set, parameters::parameter_values(set),
                                                      {.id = "other", .method_id = "eem"}),
                    std::invalid_argument);
}
