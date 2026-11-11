#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <chargefw/features/topology_features.h>
#include <chargefw/parameters/classification/parameter_classifier.h>
#include <chargefw/parameters/models/atom_parameters.h>
#include <chargefw/parameters/models/bond_parameters.h>
#include <chargefw/parameters/models/common_parameters.h>
#include <chargefw/parameters/models/parameter_set.h>
#include <chargefw/parameters/models/parameter_set_metadata.h>
#include <chargefw/parameters/models/parameter_view.h>

#include <snitch/snitch.hpp>
#include <type_traits>

namespace features = chargefw::features;
namespace parameters = chargefw::parameters;

namespace {

auto make_parameter_set() -> parameters::ParameterSet {
    return parameters::ParameterSet{
        {.id = "test-view", .method_id = "test", .name = "Test parameter view"},
        parameters::CommonParameters{{{.name = "common", .value = 4.0}}},
        parameters::AtomParameters{{{.key = chargefw::test::plain_atom_key(1),
                                     .parameters = {{.name = "atom", .value = 1.0}}},
                                    {.key = chargefw::test::plain_atom_key(8),
                                     .parameters = {{.name = "atom", .value = 2.0}}}}},
        parameters::BondParameters{{{.key = chargefw::test::plain_bond_key(1, 8),
                                     .parameters = {{.name = "bond", .value = 3.0}}}}}};
}

} // namespace

static_assert(!std::is_constructible_v<parameters::ParameterView, parameters::ParameterSet&&,
                                       const parameters::ParameterClassification&>);

TEST_CASE("parameter view exposes immutable classified parameter access", "[parameters][view]") {
    const auto water = chargefw::test::make_water_graph();
    const auto parameter_set = make_parameter_set();
    const features::TopologyFeatures topology{water};
    const auto classification = parameters::classify_parameters(topology, parameter_set);
    const parameters::ParameterView view{parameter_set, classification};

    CHECK(view.parameter_set().id() == "test-view");
    const auto atom_indices = view.classification().atom().parameter_entry_indices();
    REQUIRE(atom_indices.size() == 3);
    CHECK(atom_indices[0] == 1);
    CHECK(atom_indices[1] == 0);
    CHECK(atom_indices[2] == 0);
    CHECK(view.common("common") == 4.0);
    CHECK(view.atom("atom")[0] == 2.0);
    CHECK(view.atom("atom")[1] == 1.0);
    CHECK(view.atom("atom")[2] == 1.0);
    CHECK(view.bond("bond")[0] == 3.0);
    CHECK(view.bond("bond")[1] == 3.0);
}
