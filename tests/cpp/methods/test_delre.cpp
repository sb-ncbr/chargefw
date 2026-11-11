#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <chargefw/core/molecule.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/methods/calculation_input.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/methods/method_registry.h>
#include <chargefw/parameters/classification/parameter_classification.h>
#include <chargefw/parameters/classification/parameter_classifier.h>
#include <chargefw/parameters/models/atom_parameters.h>
#include <chargefw/parameters/models/bond_parameters.h>
#include <chargefw/parameters/models/parameter_key.h>
#include <chargefw/parameters/models/parameter_set.h>
#include <chargefw/parameters/models/parameter_set_metadata.h>
#include <chargefw/parameters/models/parameter_view.h>

#include <cmath>
#include <cstddef>
#include <snitch/snitch.hpp>
#include <vector>

namespace core = chargefw::core;
namespace features = chargefw::features;
namespace methods = chargefw::methods;
namespace parameters = chargefw::parameters;

namespace {

auto make_wildcard_parameter_set() -> parameters::ParameterSet {
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{
            .id = "wildcard-delre", .method_id = "delre", .name = "Wildcard DelRe parameters"},
        {},
        parameters::AtomParameters{{{.key = chargefw::test::plain_atom_key(8),
                                     .parameters = {{.name = "delta", .value = 2.0}}},
                                    {.key = chargefw::test::plain_atom_key(1),
                                     .parameters = {{.name = "delta", .value = 1.0}}}}},
        parameters::BondParameters{{{.key = chargefw::test::plain_bond_key(0, 1),
                                     .parameters = {{.name = "eps", .value = 1.0},
                                                    {.name = "gammaA", .value = 0.2},
                                                    {.name = "gammaB", .value = 0.1}}}}}};
}

auto calculate_delre(const methods::Method& method, const core::Molecule& molecule,
                     const parameters::ParameterSet& parameter_set,
                     const parameters::ParameterClassification& classification)
    -> chargefw::charges::AtomicCharges {
    const features::PreparedMolecule prepared{molecule};
    const auto method_options = methods::make_default_options(method.option_schema());
    const parameters::ParameterView parameter_view{parameter_set, classification};
    const methods::CalculationInput input{
        prepared, method_options, core::total_formal_charge(molecule), nullptr, &parameter_view};

    return method.calculate(input);
}

} // namespace

TEST_CASE("DelRe wildcard bond keys preserve physical endpoint orientation", "[methods][delre]") {
    const auto* delre = methods::method_registry().find("delre");
    REQUIRE(delre != nullptr);

    const auto parameter_set = make_wildcard_parameter_set();
    const auto molecule = chargefw::test::make_water_graph();
    const auto reversed_molecule = chargefw::test::flip_bond_directions(molecule);
    const auto classification =
        parameters::classify_parameters(features::TopologyFeatures{molecule}, parameter_set);
    const auto reversed_classification = parameters::classify_parameters(
        features::TopologyFeatures{reversed_molecule}, parameter_set);

    const auto charges = calculate_delre(*delre, molecule, parameter_set, classification);
    const auto reversed_charges =
        calculate_delre(*delre, reversed_molecule, parameter_set, reversed_classification);

    REQUIRE(charges.size() == reversed_charges.size());
    for (std::size_t atom_index = 0; atom_index < charges.size(); ++atom_index) {
        CHECK(std::abs(charges[atom_index] - reversed_charges[atom_index]) < 1.0e-12);
    }
}
