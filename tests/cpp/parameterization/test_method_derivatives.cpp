#include "parameterization/method_derivatives.h"
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
#include <chargefw/methods/method_options.h>
#include <chargefw/methods/method_registry.h>
#include <chargefw/parameterization/normal_equations.h>
#include <chargefw/parameterization/parameter_vector.h>
#include <chargefw/parameters/classification/parameter_classifier.h>
#include <chargefw/parameters/models/parameter_set.h>
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

// Exact derivatives of the methods that have them (parameterization/method_derivatives.h). The
// molecules use plain H and O types, so atom entry 0 is hydrogen, entry 1 oxygen, and one bond
// entry covers every bond.

namespace derivatives = chargefw::parameterization::detail;
namespace methods = chargefw::methods;
namespace parameterization = chargefw::parameterization;
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
    const auto classification = parameters::classify_parameters(prepared.topology(), parameter_set);
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
    -> std::unique_ptr<derivatives::LinearizedCharges> {
    auto result =
        with_input(molecule, method_id, parameter_set, sources, references,
                   [](const methods::Method& selected, const methods::CalculationInput& input) {
                       const auto* found = derivatives::derivatives_of(selected);
                       REQUIRE(found != nullptr);
                       return found->linearize(input);
                   });
    REQUIRE(result != nullptr);
    return result;
}

// Columns with one per slot of the parameter set's layout: atom entry 0 is hydrogen, entry 1
// oxygen, and one bond entry covers every bond.
auto slot_columns(const core::Molecule& molecule, const methods::Method& selected,
                  const parameters::ParameterSet& parameter_set) -> derivatives::ParameterColumns {
    const auto layout = parameterization::parameter_layout(parameter_set);
    const auto slot_of = [&layout](const parameterization::ParameterScope scope,
                                   const std::size_t entry, const std::string_view name) {
        const auto found =
            std::ranges::find_if(layout, [&](const parameterization::ParameterSlot& slot) {
                return slot.scope == scope && slot.entry_index == entry && slot.name == name;
            });
        REQUIRE(found != layout.end());
        return static_cast<std::size_t>(found - layout.begin());
    };
    const auto& requirements = selected.requirements();
    auto columns =
        derivatives::ParameterColumns{.count = layout.size(), .atom = {}, .bond = {}, .common = {}};
    for (const auto name : requirements.atom_parameters) {
        auto& assigned = columns.atom.emplace_back();
        for (std::size_t atom = 0; atom < molecule.atom_count(); ++atom) {
            const auto entry =
                molecule.atom(atom).atomic_number() == 1 ? std::size_t{0} : std::size_t{1};
            assigned.push_back(slot_of(parameterization::ParameterScope::atom, entry, name));
        }
    }
    for (const auto name : requirements.bond_parameters) {
        columns.bond.emplace_back(molecule.bond_count(),
                                  slot_of(parameterization::ParameterScope::bond, 0, name));
    }
    for (const auto name : requirements.common_parameters) {
        columns.common.push_back(slot_of(parameterization::ParameterScope::common, 0, name));
    }
    return columns;
}

// Column `column` of a jacobian() result for `atom_count` atoms.
auto column_of(const std::vector<double>& jacobian, const std::size_t atom_count,
               const std::size_t column) -> std::span<const double> {
    return std::span{jacobian}.subspan(column * atom_count, atom_count);
}

auto shifted(const parameters::ParameterSet& parameter_set, const std::size_t slot,
             const double amount) -> parameters::ParameterSet {
    auto values = parameterization::parameter_values(parameter_set);
    values[slot] += amount;
    return parameterization::with_parameter_values(parameter_set, values, {.id = "shifted"});
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

TEST_CASE("Exact derivatives belong to the methods that have them", "[methods][derivatives]") {
    for (const auto& [method_id, parameter_set] : cases()) {
        CAPTURE(method_id);
        CHECK(chargefw::parameterization::has_exact_derivatives(method(method_id)));
    }
    // PEOE iterates instead of solving a system; it has no exact derivatives.
    CHECK_FALSE(chargefw::parameterization::has_exact_derivatives(method("peoe")));
}

TEST_CASE("Jacobian matches central differences of every slot", "[methods][sqe][derivatives]") {
    const auto molecule = make_water_and_hydroxide();
    constexpr auto step = 1e-6;
    for (const auto& [method_id, parameter_set] : cases()) {
        const auto& selected = method(method_id);
        REQUIRE(chargefw::parameterization::has_exact_derivatives(selected));
        for (const auto sources : source_choices) {
            const auto linearized = linearize(molecule, method_id, parameter_set, sources);
            const auto columns = slot_columns(molecule, selected, parameter_set);
            const auto jacobian = linearized->jacobian(columns);
            REQUIRE(jacobian.size() == molecule.atom_count() * columns.count);

            const auto layout = parameterization::parameter_layout(parameter_set);
            for (std::size_t slot = 0; slot < layout.size(); ++slot) {
                CAPTURE(method_id, sources.size(), layout[slot].name, layout[slot].entry_index);
                const auto exact = column_of(jacobian, molecule.atom_count(), slot);
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

// Columns for `selected` on `molecule`: values share columns across parameters and atoms, the
// second atom's values belong to none, and the last column has no values. With `separate`, every
// value has a column of its own instead.
auto grouped_columns(const core::Molecule& molecule, const methods::Method& selected,
                     const bool separate) -> derivatives::ParameterColumns {
    const auto& requirements = selected.requirements();
    auto columns = derivatives::ParameterColumns{.count = 4, .atom = {}, .bond = {}, .common = {}};
    auto next = std::size_t{0};
    const auto column = [&next, separate](const std::size_t grouped) {
        return separate ? next++ : grouped;
    };
    for (std::size_t channel = 0; channel < requirements.atom_parameters.size(); ++channel) {
        auto& assigned = columns.atom.emplace_back();
        for (std::size_t atom = 0; atom < molecule.atom_count(); ++atom) {
            assigned.push_back(atom == 1 && !separate ? derivatives::ParameterColumns::none
                                                      : column((channel + atom) % 3));
        }
    }
    for (std::size_t channel = 0; channel < requirements.bond_parameters.size(); ++channel) {
        auto& assigned = columns.bond.emplace_back();
        for (std::size_t bond = 0; bond < molecule.bond_count(); ++bond) {
            assigned.push_back(column(bond % 2));
        }
    }
    for (std::size_t channel = 0; channel < requirements.common_parameters.size(); ++channel) {
        columns.common.push_back(column(2));
    }
    if (separate) {
        columns.count = next;
    }
    return columns;
}

// Checks that each grouped column is the sum of the separate columns of the values it groups.
auto check_grouping(const derivatives::LinearizedCharges& linearized,
                    const core::Molecule& molecule, const methods::Method& selected) -> void {
    const auto atom_count = molecule.atom_count();
    const auto grouped = grouped_columns(molecule, selected, false);
    const auto separate = grouped_columns(molecule, selected, true);
    const auto actual = linearized.jacobian(grouped);
    const auto values = linearized.jacobian(separate);
    REQUIRE(actual.size() == atom_count * grouped.count);
    REQUIRE(values.size() == atom_count * separate.count);

    auto expected = std::vector<double>(actual.size(), 0.0);
    const auto add = [&](const std::size_t group, const std::size_t value) {
        if (group != derivatives::ParameterColumns::none) {
            for (std::size_t atom = 0; atom < atom_count; ++atom) {
                expected[group * atom_count + atom] += values[value * atom_count + atom];
            }
        }
    };
    for (std::size_t channel = 0; channel < grouped.atom.size(); ++channel) {
        for (std::size_t atom = 0; atom < atom_count; ++atom) {
            add(grouped.atom[channel][atom], separate.atom[channel][atom]);
        }
    }
    for (std::size_t channel = 0; channel < grouped.bond.size(); ++channel) {
        for (std::size_t bond = 0; bond < grouped.bond[channel].size(); ++bond) {
            add(grouped.bond[channel][bond], separate.bond[channel][bond]);
        }
    }
    for (std::size_t channel = 0; channel < grouped.common.size(); ++channel) {
        add(grouped.common[channel], separate.common[channel]);
    }
    for (std::size_t index = 0; index < actual.size(); ++index) {
        CAPTURE(index);
        CHECK(std::abs(actual[index] - expected[index]) <=
              1e-12 * (1.0 + std::abs(expected[index])));
    }
}

TEST_CASE("Jacobian columns add up the values they group", "[methods][derivatives]") {
    const auto molecule = make_water_and_hydroxide();
    for (const auto& [method_id, parameter_set] : cases()) {
        const auto& selected = method(method_id);
        for (const auto sources : source_choices) {
            CAPTURE(method_id, sources.size());
            check_grouping(*linearize(molecule, method_id, parameter_set, sources), molecule,
                           selected);
        }
        // Without bonds: SQE's bond parameter has no values, and SQE+qp still normalizes q0.
        const auto separate_atoms = make_separate_atoms();
        check_grouping(*linearize(separate_atoms, method_id, parameter_set, {}), separate_atoms,
                       selected);
    }
}

TEST_CASE("SQE+qp q0 changes keep every bonded part's charge, and references bypass q0",
          "[methods][sqe][derivatives]") {
    const auto molecule = make_water_and_hydroxide();
    const auto parameter_set = cases()[2].parameter_set;
    const auto& selected = method("sqeqp");
    const auto columns = slot_columns(molecule, selected, parameter_set);
    const auto q0_slots = std::array{columns.atom[3][1], columns.atom[3][0]}; // H, then O.

    // The total charge of each bonded part does not depend on any q0, while single charges do.
    const auto jacobian = linearize(molecule, "sqeqp", parameter_set, {})->jacobian(columns);
    for (const auto slot : q0_slots) {
        const auto change = column_of(jacobian, molecule.atom_count(), slot);
        CHECK(std::abs(change[0] + change[1] + change[2]) < 1e-12);
        CHECK(std::abs(change[3] + change[4]) < 1e-12);
        CHECK(std::abs(change[0]) > 1e-3);
    }

    const auto references = std::vector{-0.6, 0.3, 0.3, -1.1, 0.1};
    const auto bypassed =
        linearize(molecule, "sqeqp", parameter_set, {}, references)->jacobian(columns);
    for (const auto slot : q0_slots) {
        CHECK(all_zero(column_of(bypassed, molecule.atom_count(), slot)));
    }
}

TEST_CASE("SQE-family linearization covers molecules without bonds",
          "[methods][sqe][derivatives]") {
    // Without bonds, SQE charges stay zero; SQE+qp keeps each atom's normalized q0, which no q0
    // change can move.
    const auto separate = make_separate_atoms();
    for (const auto& [method_id, parameter_set] : cases()) {
        const auto& selected = method(method_id);
        if (!selected.requirements().bonds) {
            continue;
        }
        CAPTURE(method_id);
        const auto linearized = linearize(separate, method_id, parameter_set, {});
        CHECK(all_zero(linearized->jacobian(slot_columns(separate, selected, parameter_set))));
    }
}

TEST_CASE("Linearization outlives its input and rejects malformed columns",
          "[methods][derivatives]") {
    const auto all_cases = cases();
    const auto molecule = make_water_and_hydroxide();
    for (std::size_t index = 0; index < all_cases.size(); ++index) {
        const auto& [method_id, parameter_set] = all_cases[index];
        const auto& selected = method(method_id);
        CAPTURE(method_id);
        const auto charges = charges_of(molecule, method_id, parameter_set, source);
        // The molecule, the parameter set, and the input objects are destroyed before any use.
        const auto linearized = [&method_id, index] {
            const auto temporary_molecule = make_water_and_hydroxide();
            const auto temporary_parameters = cases()[index].parameter_set;
            return linearize(temporary_molecule, method_id, temporary_parameters, source);
        }();
        CHECK(std::ranges::equal(linearized->charges(), charges));
        const auto columns = slot_columns(molecule, selected, parameter_set);
        const auto jacobian = linearized->jacobian(columns);
        for (std::size_t slot = 0; slot < columns.count; ++slot) {
            CAPTURE(slot);
            CHECK_FALSE(all_zero(column_of(jacobian, molecule.atom_count(), slot)));
        }

        for (const auto channel : {std::size_t{0}, columns.atom.size() - 1}) {
            CAPTURE(channel);
            auto outside = columns;
            outside.atom[channel][0] = outside.count;
            CHECK_THROWS_AS(linearized->jacobian(outside), std::invalid_argument);
            auto short_channel = columns;
            short_channel.atom[channel].pop_back();
            CHECK_THROWS_AS(linearized->jacobian(short_channel), std::invalid_argument);
        }
        auto missing = columns;
        missing.atom.pop_back();
        CHECK_THROWS_AS(linearized->jacobian(missing), std::invalid_argument);
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
