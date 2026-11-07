#include "support/test_parameters.h"

#include <chargefw/charges/atomic_charges.h>
#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/position.h>
#include <chargefw/features/conformer_features.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/methods/calculation_input.h>
#include <chargefw/methods/method.h>
#include <chargefw/methods/method_derivatives.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/methods/method_registry.h>
#include <chargefw/parameters/classification/parameter_classifier.h>
#include <chargefw/parameters/models/parameter_set.h>
#include <chargefw/parameters/models/parameter_values.h>
#include <chargefw/parameters/models/parameter_view.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <memory>
#include <snitch/snitch.hpp>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

// Exact derivatives of the methods that have them (Method::linearize). The molecules use plain H
// and O types, so atom entry 0 is hydrogen, entry 1 oxygen, and one bond entry covers every bond.

namespace methods = chargefw::methods;
namespace parameters = chargefw::parameters;
namespace core = chargefw::core;

namespace {

// Water and a hydroxide ion 3 Å away: two bonded parts, one of them charged.
auto make_water_and_hydroxide() -> core::Molecule {
    std::vector atoms{core::Atom{8, 0, "O1"}, core::Atom{1, 0, "H1"}, core::Atom{1, 0, "H2"},
                      core::Atom{8, -1, "O2"}, core::Atom{1, 0, "H3"}};
    std::vector bonds{core::Bond{0, 1}, core::Bond{0, 2}, core::Bond{3, 4}};
    std::vector positions{core::Position{.x = 0.0, .y = 0.0, .z = 0.0},
                          core::Position{.x = 0.96, .y = 0.0, .z = 0.0},
                          core::Position{.x = -0.24, .y = 0.93, .z = 0.0},
                          core::Position{.x = 3.0, .y = 0.2, .z = 0.1},
                          core::Position{.x = 3.97, .y = 0.3, .z = 0.1}};
    std::vector conformers{core::Conformer{std::move(positions), "model-1"}};
    return core::Molecule{std::move(atoms), std::move(bonds), std::move(conformers), "two-parts"};
}

auto make_separate_atoms() -> core::Molecule {
    std::vector atoms{core::Atom{8, 0, "O"}, core::Atom{1, 0, "H"}};
    std::vector positions{core::Position{.x = 0.0, .y = 0.0, .z = 0.0},
                          core::Position{.x = 2.0, .y = 0.0, .z = 0.0}};
    std::vector conformers{core::Conformer{std::move(positions), "model-1"}};
    return core::Molecule{std::move(atoms), {}, std::move(conformers), "no-bonds"};
}

struct Case {
    std::string method_id;
    parameters::ParameterSet parameter_set;
};

auto cases() -> std::vector<Case> {
    return {
        {"sqe", chargefw::test::make_sqe_ho_parameters("sqe", {.hydrogen_width = 0.8})},
        {"sqeq0", chargefw::test::make_sqe_ho_parameters("sqeq0", {.hydrogen_width = 0.8})},
        {"sqeqp", chargefw::test::make_sqe_ho_parameters(
                      "sqeqp", {.hydrogen_width = 0.8, .hydrogen_q0 = 0.1, .oxygen_q0 = -0.3})},
        {"eem", chargefw::test::make_eem_ho_parameters(0.5)},
    };
}

auto method(const std::string_view method_id) -> const methods::Method& {
    const auto* found = methods::method_registry().find(method_id);
    if (found == nullptr) {
        throw std::logic_error{"method is not registered"};
    }
    return *found;
}

// Calls `operation` with the input of `molecule` and `parameter_set`; every object the input
// borrows lives only during the call.
template <typename Operation>
auto with_input(const core::Molecule& molecule, const std::string_view method_id,
                const parameters::ParameterSet& parameter_set,
                const std::span<const methods::FixedPointSource> sources,
                const std::span<const double> references, Operation operation) {
    const chargefw::features::PreparedMolecule prepared{molecule};
    const chargefw::features::ConformerFeatures geometry{molecule, 0};
    const auto classification =
        parameters::classify_parameters(molecule, prepared.topology(), parameter_set);
    const parameters::ParameterView view{parameter_set, classification};
    const auto options = methods::make_default_options(method(method_id).option_schema());
    const methods::CalculationInput input{prepared, options, -1.0,      &geometry,
                                          &view,    sources, references};
    return operation(method(method_id), input);
}

auto charges_of(const core::Molecule& molecule, const std::string_view method_id,
                const parameters::ParameterSet& parameter_set,
                const std::span<const methods::FixedPointSource> sources) -> std::vector<double> {
    return with_input(molecule, method_id, parameter_set, sources, {},
                      [](const methods::Method& selected, const methods::CalculationInput& input) {
                          const auto charges = selected.calculate(input);
                          return std::vector<double>(charges.values().begin(),
                                                     charges.values().end());
                      });
}

auto linearize(const core::Molecule& molecule, const std::string_view method_id,
               const parameters::ParameterSet& parameter_set,
               const std::span<const methods::FixedPointSource> sources,
               const std::span<const double> references = {})
    -> std::unique_ptr<methods::LinearizedCharges> {
    auto result =
        with_input(molecule, method_id, parameter_set, sources, references,
                   [](const methods::Method& selected, const methods::CalculationInput& input) {
                       return selected.linearize(input);
                   });
    REQUIRE(result != nullptr);
    return result;
}

// The derivative of every charge with respect to one slot of the parameter set: the sum of the
// charge's gradient over every atom or bond the slot's entry matches.
auto slot_derivatives(const core::Molecule& molecule, const methods::Method& selected,
                      const std::span<const methods::ParameterGradient> charge_gradients,
                      const parameters::ParameterSlot& slot) -> std::vector<double> {
    const auto& requirements = selected.requirements();
    const auto channel = [&slot](const std::vector<std::string_view>& names) {
        return static_cast<std::size_t>(std::ranges::find(names, slot.name) - names.begin());
    };
    auto result = std::vector<double>{};
    for (const auto& gradient : charge_gradients) {
        auto sum = 0.0;
        if (slot.scope == parameters::ParameterScope::atom) {
            const auto element = slot.entry_index == 0 ? 1 : 8;
            for (std::size_t atom = 0; atom < molecule.atom_count(); ++atom) {
                if (molecule.atom(atom).atomic_number() == element) {
                    sum += gradient.atom[channel(requirements.atom_parameters)][atom];
                }
            }
        } else if (slot.scope == parameters::ParameterScope::bond) {
            for (const auto value : gradient.bond[channel(requirements.bond_parameters)]) {
                sum += value;
            }
        } else {
            sum = gradient.common[channel(requirements.common_parameters)];
        }
        result.push_back(sum);
    }
    return result;
}

auto charge_gradients(const methods::LinearizedCharges& linearized)
    -> std::vector<methods::ParameterGradient> {
    const auto atom_count = linearized.charges().size();
    auto result = std::vector<methods::ParameterGradient>{};
    for (std::size_t atom = 0; atom < atom_count; ++atom) {
        auto cotangent = std::vector<double>(atom_count, 0.0);
        cotangent[atom] = 1.0;
        result.push_back(linearized.gradient(cotangent));
    }
    return result;
}

auto shifted(const parameters::ParameterSet& parameter_set, const std::size_t slot,
             const double amount) -> parameters::ParameterSet {
    auto values = parameters::parameter_values(parameter_set);
    values[slot] += amount;
    return parameters::with_parameter_values(parameter_set, values, {.id = "shifted"});
}

auto all_zero(const std::span<const double> values) -> bool {
    return std::ranges::all_of(values, [](const double value) { return value == 0.0; });
}

const std::array source{methods::FixedPointSource{
    .position = core::Position{.x = 1.5, .y = -1.2, .z = 0.4}, .charge = 0.7}};
const auto source_choices = std::array{std::span<const methods::FixedPointSource>{},
                                       std::span<const methods::FixedPointSource>{source}};

} // namespace

TEST_CASE("Linearization reproduces the calculated charges", "[methods][sqe][derivatives]") {
    const auto molecule = make_water_and_hydroxide();
    for (const auto& [method_id, parameter_set] : cases()) {
        for (const auto sources : source_choices) {
            CAPTURE(method_id, sources.size());
            const auto linearized = linearize(molecule, method_id, parameter_set, sources);
            const auto charges = charges_of(molecule, method_id, parameter_set, sources);
            const auto linearized_charges = linearized->charges();
            CHECK(std::ranges::equal(linearized_charges, charges));
        }
    }
}

TEST_CASE("Gradients match central differences of every slot", "[methods][sqe][derivatives]") {
    const auto molecule = make_water_and_hydroxide();
    constexpr auto step = 1e-6;
    for (const auto& [method_id, parameter_set] : cases()) {
        const auto& selected = method(method_id);
        for (const auto sources : source_choices) {
            REQUIRE(selected.requirements().supports_exact_derivatives);
            const auto linearized = linearize(molecule, method_id, parameter_set, sources);
            const auto gradients = charge_gradients(*linearized);
            // The channels are the method's parameters, in its order.
            REQUIRE(gradients.front().atom.size() ==
                    selected.requirements().atom_parameters.size());
            REQUIRE(gradients.front().bond.size() ==
                    selected.requirements().bond_parameters.size());
            REQUIRE(gradients.front().common.size() ==
                    selected.requirements().common_parameters.size());

            const auto layout = parameters::parameter_layout(parameter_set);
            for (std::size_t slot = 0; slot < layout.size(); ++slot) {
                CAPTURE(method_id, sources.size(), layout[slot].name, layout[slot].entry_index);
                const auto exact = slot_derivatives(molecule, selected, gradients, layout[slot]);
                const auto up =
                    charges_of(molecule, method_id, shifted(parameter_set, slot, step), sources);
                const auto down =
                    charges_of(molecule, method_id, shifted(parameter_set, slot, -step), sources);
                for (std::size_t atom = 0; atom < exact.size(); ++atom) {
                    const auto central = (up[atom] - down[atom]) / (2.0 * step);
                    CHECK(std::abs(exact[atom] - central) <= 1e-7 * (1.0 + std::abs(central)));
                }
            }
        }
    }
}

TEST_CASE("SQE+qp q0 changes keep every bonded part's charge, and references bypass q0",
          "[methods][sqe][derivatives]") {
    const auto molecule = make_water_and_hydroxide();
    const auto parameter_set = cases()[2].parameter_set;
    const auto linearized = linearize(molecule, "sqeqp", parameter_set, {});
    constexpr std::size_t q0 = 3; // The channel order of SQE+qp's requirements().

    // The total charge of each bonded part does not depend on any q0, while single charges do.
    const auto water = linearized->gradient(std::vector{1.0, 1.0, 1.0, 0.0, 0.0});
    const auto hydroxide = linearized->gradient(std::vector{0.0, 0.0, 0.0, 1.0, 1.0});
    for (std::size_t atom = 0; atom < molecule.atom_count(); ++atom) {
        CHECK(std::abs(water.atom[q0][atom]) < 1e-12);
        CHECK(std::abs(hydroxide.atom[q0][atom]) < 1e-12);
    }
    const auto oxygen = linearized->gradient(std::vector{1.0, 0.0, 0.0, 0.0, 0.0});
    CHECK(std::ranges::any_of(oxygen.atom[q0],
                              [](const double value) { return std::abs(value) > 1e-3; }));

    const auto references = std::vector{-0.6, 0.3, 0.3, -1.1, 0.1};
    const auto bypassed = linearize(molecule, "sqeqp", parameter_set, {}, references);
    CHECK(all_zero(bypassed->gradient(std::vector{1.0, 0.0, -1.0, 0.5, 0.2}).atom[q0]));
}

TEST_CASE("SQE-family linearization covers molecules without bonds",
          "[methods][sqe][derivatives]") {
    // Without bonds, SQE charges stay zero; SQE+qp keeps each atom's normalized q0, which no q0
    // change can move.
    const auto separate = make_separate_atoms();
    for (const auto& [method_id, parameter_set] : cases()) {
        if (!method(method_id).requirements().bonds) {
            continue;
        }
        CAPTURE(method_id);
        const auto linearized = linearize(separate, method_id, parameter_set, {});
        const auto gradient = linearized->gradient(std::vector{1.0, -2.0});
        CHECK(gradient.atom.size() == method(method_id).requirements().atom_parameters.size());
        for (const auto& channel : gradient.atom) {
            CHECK(all_zero(channel));
        }
    }
}

TEST_CASE("Linearization outlives its input and rejects malformed cotangents",
          "[methods][derivatives]") {
    const auto all_cases = cases();
    for (std::size_t index = 0; index < all_cases.size(); ++index) {
        const auto& [method_id, parameter_set] = all_cases[index];
        CAPTURE(method_id);
        const auto charges =
            charges_of(make_water_and_hydroxide(), method_id, parameter_set, source);
        // The molecule, the parameter set, and the input objects are destroyed before any use.
        const auto linearized = [&method_id, index] {
            const auto molecule = make_water_and_hydroxide();
            const auto temporary_parameters = cases()[index].parameter_set;
            return linearize(molecule, method_id, temporary_parameters, source);
        }();
        CHECK(std::ranges::equal(linearized->charges(), charges));
        const auto gradient = linearized->gradient(std::vector{0.4, -1.1, 0.7, 0.25, -0.9});
        for (const auto& channel : gradient.atom) {
            CHECK_FALSE(all_zero(channel));
        }

        CHECK_THROWS_AS(linearized->gradient(std::vector{1.0}), std::invalid_argument);
        CHECK_THROWS_AS(linearized->gradient(std::vector{1.0, 0.0, std::nan(""), 0.0, 0.0}),
                        std::invalid_argument);
        if (method_id == "sqeqp") {
            // Summing each bonded part's q0 gradient overflows; the gradient is not finite.
            CHECK_THROWS_AS(linearized->gradient(std::vector(5, 1e308)), std::runtime_error);
        }
    }
}

TEST_CASE("Linearization fails where calculation fails", "[methods][derivatives]") {
    // Two atoms 1 Å apart with kappa = 1 and both B = 1 make the EEM system singular.
    std::vector atoms{core::Atom{1, 0, "H"}, core::Atom{8, 0, "O"}};
    std::vector positions{core::Position{.x = 0.0, .y = 0.0, .z = 0.0},
                          core::Position{.x = 1.0, .y = 0.0, .z = 0.0}};
    std::vector conformers{core::Conformer{std::move(positions), "model-1"}};
    const core::Molecule molecule{std::move(atoms), {}, std::move(conformers), "singular"};
    const parameters::ParameterSet parameter_set{
        parameters::ParameterSetMetadata{.id = "singular", .method_id = "eem", .name = "Singular"},
        parameters::CommonParameters{{{.name = "kappa", .value = 1.0}}},
        parameters::AtomParameters{
            {{.key = chargefw::test::plain_atom_key(1),
              .parameters = {{.name = "A", .value = 1.0}, {.name = "B", .value = 1.0}}},
             {.key = chargefw::test::plain_atom_key(8),
              .parameters = {{.name = "A", .value = 2.0}, {.name = "B", .value = 1.0}}}}}};

    CHECK_THROWS_AS(charges_of(molecule, "eem", parameter_set, {}), std::invalid_argument);
    CHECK_THROWS_AS(linearize(molecule, "eem", parameter_set, {}), std::invalid_argument);
}
