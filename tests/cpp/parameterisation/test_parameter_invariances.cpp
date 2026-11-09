#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <chargefw/charges/atomic_charges.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/position.h>
#include <chargefw/features/conformer_features.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/methods/calculation_input.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/methods/method_registry.h>
#include <chargefw/parameterisation/parameter_vector.h>
#include <chargefw/parameters/classification/parameter_classifier.h>
#include <chargefw/parameters/models/parameter_set.h>
#include <chargefw/parameters/models/parameter_view.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <snitch/snitch.hpp>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

// The fitting profiles fix one slot per invariance: a change of parameter values that changes no
// charge. These tests check each declared invariance on the methods themselves, with and without
// fixed point sources, so that a profile never fixes a slot the data could determine.

namespace methods = chargefw::methods;
namespace parameterisation = chargefw::parameterisation;
namespace parameters = chargefw::parameters;

namespace {

using Transform = std::function<double(const parameterisation::ParameterSlot&, double)>;

auto calculate(const chargefw::core::Molecule& molecule, const std::string_view method_id,
               const parameters::ParameterSet& parameter_set,
               const std::span<const methods::FixedPointSource> sources)
    -> chargefw::charges::AtomicCharges {
    const chargefw::features::PreparedMolecule prepared{molecule};
    const chargefw::features::ConformerFeatures geometry{molecule, 0};
    const auto classification =
        parameters::classify_parameters(molecule, prepared.topology(), parameter_set);
    const parameters::ParameterView view{parameter_set, classification};
    const auto* method = methods::method_registry().find(method_id);
    if (method == nullptr) {
        throw std::logic_error{"method is not registered"};
    }
    const auto options = methods::make_default_options(method->option_schema());
    const methods::CalculationInput input{prepared, options, 0.0, &geometry, &view, sources};
    return method->calculate(input);
}

auto transformed(const parameters::ParameterSet& parameter_set, const Transform& transform)
    -> parameters::ParameterSet {
    const auto layout = parameterisation::parameter_layout(parameter_set);
    auto values = parameterisation::parameter_values(parameter_set);
    for (std::size_t index = 0; index < layout.size(); ++index) {
        values[index] = transform(layout[index], values[index]);
    }
    return parameterisation::with_parameter_values(parameter_set, values, {.id = "transformed"});
}

auto max_difference(const chargefw::charges::AtomicCharges& first,
                    const chargefw::charges::AtomicCharges& second) -> double {
    REQUIRE(first.size() == second.size());
    auto difference = 0.0;
    for (std::size_t index = 0; index < first.size(); ++index) {
        difference = std::max(difference, std::abs(first[index] - second[index]));
    }
    return difference;
}

auto add_to(const std::string_view name, const double amount) -> Transform {
    return [name, amount](const parameterisation::ParameterSlot& slot, const double value) {
        return slot.name == name ? value + amount : value;
    };
}

// Changes only the oxygen entry, the second atom entry of the H/O test sets.
auto add_to_oxygen(const std::string_view name, const double amount) -> Transform {
    return [name, amount](const parameterisation::ParameterSlot& slot, const double value) {
        return slot.name == name && slot.entry_index == 1 ? value + amount : value;
    };
}

struct InvarianceCase {
    std::string label;
    std::string method_id;
    parameters::ParameterSet parameter_set;
    Transform invariant;
    Transform contrast;
};

auto cases() -> std::vector<InvarianceCase> {
    const auto eem = chargefw::test::make_eem_ho_parameters(0.5);
    const auto sqe = chargefw::test::make_sqe_ho_parameters("sqe");
    const auto sqeq0 = chargefw::test::make_sqe_ho_parameters("sqeq0");
    const auto sqeqp =
        chargefw::test::make_sqe_ho_parameters("sqeqp", {.hydrogen_q0 = 0.1, .oxygen_q0 = -0.3});
    const auto scale_all = [](const parameterisation::ParameterSlot& /*slot*/, const double value) {
        return 1.7 * value;
    };
    const auto flip_hydrogen_width = [](const parameterisation::ParameterSlot& slot,
                                        const double value) {
        return slot.name == "width" && slot.entry_index == 0 ? -value : value;
    };
    return {
        {"EEM: same constant added to every A", "eem", eem, add_to("A", 0.7),
         add_to_oxygen("A", 0.7)},
        {"EEM: A, B, and kappa scaled together", "eem", eem, scale_all, add_to_oxygen("B", 0.7)},
        {"SQE: same constant added to every electronegativity", "sqe", sqe,
         add_to("electronegativity", 0.7), add_to_oxygen("electronegativity", 0.7)},
        {"SQE+q0: same constant added to every electronegativity", "sqeq0", sqeq0,
         add_to("electronegativity", 0.7), add_to_oxygen("electronegativity", 0.7)},
        {"SQE+qp: same constant added to every electronegativity", "sqeqp", sqeqp,
         add_to("electronegativity", 0.7), add_to_oxygen("electronegativity", 0.7)},
        {"SQE+qp: same constant added to every q0", "sqeqp", sqeqp, add_to("q0", 0.7),
         add_to_oxygen("q0", 0.7)},
        {"SQE: sign of a width", "sqe", sqe, flip_hydrogen_width, add_to_oxygen("width", 0.7)},
        {"SQE+qp: sign of a width", "sqeqp", sqeqp, flip_hydrogen_width,
         add_to_oxygen("width", 0.7)},
    };
}

} // namespace

TEST_CASE("declared parameter invariances leave every charge unchanged",
          "[methods][parameter-invariances]") {
    const auto water = chargefw::test::make_water();
    const std::array source{methods::FixedPointSource{
        .position = chargefw::core::Position{.x = 2.0, .y = 0.5, .z = 0.3}, .charge = 0.5}};
    const auto source_choices = std::array{std::span<const methods::FixedPointSource>{},
                                           std::span<const methods::FixedPointSource>{source}};

    for (const auto& [label, method_id, parameter_set, invariant, contrast] : cases()) {
        for (const auto sources : source_choices) {
            CAPTURE(label, sources.size());
            const auto original = calculate(water, method_id, parameter_set, sources);
            const auto unchanged =
                calculate(water, method_id, transformed(parameter_set, invariant), sources);
            CHECK(max_difference(original, unchanged) < 1e-12);

            // A change of one type only does change the charges, so the comparison can fail.
            const auto changed =
                calculate(water, method_id, transformed(parameter_set, contrast), sources);
            CHECK(max_difference(original, changed) > 1e-6);
        }
    }
}
