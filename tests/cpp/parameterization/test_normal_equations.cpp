#include "support/test_parameters.h"

#include <chargefw/calculation/assessment.h>
#include <chargefw/calculation/calculation.h>
#include <chargefw/calculation/execution_policy.h>
#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/core/position.h>
#include <chargefw/parameterization/normal_equations.h>
#include <chargefw/parameterization/parameter_vector.h>
#include <chargefw/parameters/models/atom_parameters.h>
#include <chargefw/parameters/models/parameter_set.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <ranges>
#include <snitch/snitch.hpp>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Normal equations from exact derivatives (parameterization::normal_equations), checked against a
// Jacobian from central differences of ordinary calculation with candidate values.

namespace calculation = chargefw::calculation;
namespace core = chargefw::core;
namespace parameterization = chargefw::parameterization;
namespace parameters = chargefw::parameters;

namespace {

auto water(const std::string& name, const std::vector<std::vector<core::Position>>& conformers)
    -> core::Molecule {
    auto models = std::vector<core::Conformer>{};
    for (std::size_t index = 0; index < conformers.size(); ++index) {
        models.emplace_back(conformers[index], "model-" + std::to_string(index + 1));
    }
    return core::Molecule{{core::Atom{8, 0, "O"}, core::Atom{1, 0, "H1"}, core::Atom{1, 0, "H2"}},
                          {core::Bond{0, 1}, core::Bond{0, 2}},
                          std::move(models),
                          name};
}

// Water with two conformers, and a hydroxide ion.
auto sqe_molecules() -> core::MoleculeCollection {
    const auto hydroxide = core::Molecule{
        {core::Atom{8, -1, "O"}, core::Atom{1, 0, "H"}},
        {core::Bond{0, 1}},
        {core::Conformer{{{.x = 0.0, .y = 0.0, .z = 0.0}, {.x = 0.97, .y = 0.1, .z = 0.0}}, "m"}},
        "hydroxide"};
    return core::MoleculeCollection{
        std::vector{water("water", {{{.x = 0.0, .y = 0.0, .z = 0.0},
                                     {.x = 0.96, .y = 0.0, .z = 0.0},
                                     {.x = -0.24, .y = 0.93, .z = 0.0}},
                                    {{.x = 0.0, .y = 0.0, .z = 0.0},
                                     {.x = 1.02, .y = 0.05, .z = 0.0},
                                     {.x = -0.31, .y = 0.90, .z = 0.1}}}),
                    hydroxide}};
}

// SQE+qp parameters stored in another order than the method declares them (q0, width, hardness,
// electronegativity), with a carbon entry that no molecule uses.
auto unusual_sqeqp_parameters() -> parameters::ParameterSet {
    const auto base = chargefw::test::make_sqe_ho_parameters(
        "sqeqp", {.hydrogen_width = 0.8, .hydrogen_q0 = 0.1, .oxygen_q0 = -0.3, .kappa = 1.2});
    auto entries = std::vector<parameters::AtomParameterEntry>{};
    for (const auto& entry : base.atom().entries()) {
        auto reversed = entry;
        std::ranges::reverse(reversed.parameters);
        entries.push_back(std::move(reversed));
    }
    entries.push_back({.key = chargefw::test::plain_atom_key(6),
                       .parameters = {{.name = "q0", .value = 0.0},
                                      {.name = "width", .value = 1.0},
                                      {.name = "hardness", .value = 10.0},
                                      {.name = "electronegativity", .value = 6.0}}});
    return parameters::ParameterSet{{.id = "unusual", .method_id = "sqeqp", .name = "Unusual"},
                                    base.common(),
                                    parameters::AtomParameters{std::move(entries)},
                                    base.bond()};
}

auto assess(core::MoleculeCollection molecules, parameters::ParameterSet parameter_set,
            std::optional<calculation::FixedIons> fixed_ions = std::nullopt,
            const calculation::ExecutionSelection selection = calculation::ExecutionSelection{
                calculation::ExecutionSelectionKind::full}) -> calculation::AssessmentResult {
    const auto method_id = std::string{parameter_set.method_id()};
    const auto parameter_set_id = std::string{parameter_set.id()};
    return calculation::assess(
        calculation::AssessmentRequest{.molecules = std::move(molecules),
                                       .parameter_sets = {std::move(parameter_set)},
                                       .method_id = method_id,
                                       .parameter_set_id = parameter_set_id,
                                       .execution_selection = selection,
                                       .fixed_ions = std::move(fixed_ions)});
}

// Every charge of every assignment, concatenated: the rows.
auto charges_of(const calculation::ExecutionPlan& plan, const parameters::ParameterSet& candidate)
    -> std::vector<double> {
    const auto result = calculation::calculate(plan, candidate);
    REQUIRE(result.calculated());
    auto rows = std::vector<double>{};
    for (const auto& assignment : result.charges->assignments()) {
        const auto values = assignment.charges.values();
        rows.insert(rows.end(), values.begin(), values.end());
    }
    return rows;
}

// The Jacobian of the rows by central differences, rows x slots, row-major.
auto differences(const calculation::ExecutionPlan& plan, const parameters::ParameterSet& candidate)
    -> std::vector<double> {
    constexpr auto step = 1e-6;
    const auto values = parameterization::parameter_values(candidate);
    const auto rows = charges_of(plan, candidate).size();
    auto jacobian = std::vector<double>(rows * values.size());
    for (std::size_t slot = 0; slot < values.size(); ++slot) {
        auto up = values;
        auto down = values;
        up[slot] += step;
        down[slot] -= step;
        const auto high =
            charges_of(plan, parameterization::with_parameter_values(candidate, up, {.id = "up"}));
        const auto low = charges_of(
            plan, parameterization::with_parameter_values(candidate, down, {.id = "down"}));
        for (std::size_t row = 0; row < rows; ++row) {
            jacobian[row * values.size() + slot] = (high[row] - low[row]) / (2.0 * step);
        }
    }
    return jacobian;
}

auto close(const std::vector<double>& actual, const std::vector<double>& expected) -> bool {
    auto scale = 1.0;
    for (const auto value : expected) {
        scale = std::max(scale, std::abs(value));
    }
    for (std::size_t index = 0; index < expected.size(); ++index) {
        if (std::abs(actual[index] - expected[index]) > 1e-6 * scale) {
            return false;
        }
    }
    return actual.size() == expected.size();
}

// Compares normal_equations() with J^T W^2 J, J^T W^2 r, and sum (w r)^2 from central differences.
auto check_against_differences(const calculation::AssessmentResult& assessment,
                               const std::vector<double>& references,
                               const std::vector<double>& weights)
    -> parameterization::NormalEquations {
    const auto& plan = *assessment.default_plan();
    const auto candidate = parameterization::with_parameter_values(
        *plan.parameter_set(), parameterization::parameter_values(*plan.parameter_set()),
        {.id = "candidate"});
    const auto result = parameterization::normal_equations(plan, candidate, references, weights);

    const auto charges = charges_of(plan, candidate);
    CHECK(result.charges == charges);

    const auto jacobian = differences(plan, candidate);
    const auto slots = parameterization::parameter_layout(candidate).size();
    auto gauss_newton = std::vector<double>(slots * slots, 0.0);
    auto gradient = std::vector<double>(slots, 0.0);
    auto sum_of_squares = 0.0;
    for (std::size_t row = 0; row < charges.size(); ++row) {
        const auto weight = weights[row] * weights[row];
        const auto residual = charges[row] - references[row];
        sum_of_squares += weight * residual * residual;
        for (std::size_t first = 0; first < slots; ++first) {
            gradient[first] += weight * jacobian[row * slots + first] * residual;
            for (std::size_t second = 0; second < slots; ++second) {
                gauss_newton[first * slots + second] +=
                    weight * jacobian[row * slots + first] * jacobian[row * slots + second];
            }
        }
    }
    CHECK(close(result.gauss_newton, gauss_newton));
    CHECK(close(result.gradient, gradient));
    CHECK(std::abs(result.sum_of_squares - sum_of_squares) <= 1e-12 * (1.0 + sum_of_squares));
    return result;
}

auto row_values(const std::size_t rows, const double start, const double step)
    -> std::vector<double> {
    auto values = std::vector<double>(rows);
    for (std::size_t row = 0; row < rows; ++row) {
        values[row] = start + step * static_cast<double>(row);
    }
    return values;
}

} // namespace

TEST_CASE("Normal equations match differences of the charges, by parameter name",
          "[calculation][derivatives]") {
    // Water with two conformers and a hydroxide ion give 3 + 3 + 2 rows.
    const auto assessment = assess(sqe_molecules(), unusual_sqeqp_parameters());
    REQUIRE(assessment.default_plan() != nullptr);
    auto weights = row_values(8, 0.5, 0.2);
    weights[4] = 0.0; // a row that is not used
    const auto result = check_against_differences(assessment, row_values(8, -0.4, 0.1), weights);

    // The carbon entry's slots have no derivatives.
    const auto layout =
        parameterization::parameter_layout(*assessment.default_plan()->parameter_set());
    const auto slots = layout.size();
    for (std::size_t slot = 0; slot < slots; ++slot) {
        if (layout[slot].scope == parameterization::ParameterScope::atom &&
            layout[slot].entry_index == 2) {
            CAPTURE(layout[slot].name);
            CHECK(result.gradient[slot] == 0.0);
            for (std::size_t other = 0; other < slots; ++other) {
                CHECK(result.gauss_newton[slot * slots + other] == 0.0);
            }
        }
    }
}

TEST_CASE("Normal equations count fixed ions as rows without derivatives",
          "[calculation][derivatives]") {
    // Water with a magnesium ion as its second atom, with a fixed charge: the water atoms are rows
    // 0, 2, and 3, the ion row 1.
    const auto molecule = core::Molecule{{core::Atom{8, 0, "O"}, core::Atom{12, 2, "Mg"},
                                          core::Atom{1, 0, "H1"}, core::Atom{1, 0, "H2"}},
                                         {core::Bond{0, 2}, core::Bond{0, 3}},
                                         {core::Conformer{{{.x = 0.0, .y = 0.0, .z = 0.0},
                                                           {.x = 3.0, .y = 0.5, .z = 0.0},
                                                           {.x = 0.96, .y = 0.0, .z = 0.0},
                                                           {.x = -0.24, .y = 0.93, .z = 0.0}},
                                                          "model-1"}},
                                         "water-with-ion"};
    const auto assessment = assess(
        core::MoleculeCollection{std::vector{molecule}},
        chargefw::test::make_eem_ho_parameters(0.5),
        calculation::FixedIons{.sources = {{.molecule_index = 0, .atom_index = 1, .charge = 2.0}}});
    REQUIRE(assessment.default_plan() != nullptr);
    const auto references = std::vector{-0.5, 1.6, 0.2, 0.3};
    const auto result = check_against_differences(assessment, references, {1.0, 3.0, 2.0, 0.5});
    CHECK(result.charges[1] == 2.0);
    // The ion's residual, (3 * (2.0 - 1.6))^2 = 1.44, is part of the sum of squares.
    CHECK(result.sum_of_squares >= 1.44);
}

TEST_CASE("Normal equations do not depend on the number of threads", "[calculation][derivatives]") {
    // More molecules than one batch of targets.
    auto molecules = std::vector<core::Molecule>{};
    for (std::size_t index = 0; index < 300; ++index) {
        const auto shift = 0.01 * static_cast<double>(index);
        molecules.push_back(
            water("water-" + std::to_string(index), {{{.x = 0.0, .y = 0.0, .z = 0.0},
                                                      {.x = 0.96 + shift, .y = 0.0, .z = 0.0},
                                                      {.x = -0.24, .y = 0.93 - shift, .z = 0.0}}}));
    }
    const auto assessment =
        assess(core::MoleculeCollection{std::move(molecules)}, unusual_sqeqp_parameters());
    const auto& plan = *assessment.default_plan();
    const auto references = row_values(900, -0.3, 0.001);
    const auto weights = row_values(900, 0.2, 0.002);
    const auto candidate = parameterization::with_parameter_values(
        *plan.parameter_set(), parameterization::parameter_values(*plan.parameter_set()),
        {.id = "candidate"});
    const auto one = parameterization::normal_equations(plan, candidate, references, weights, 1);
    const auto four = parameterization::normal_equations(plan, candidate, references, weights, 4);
    CHECK(one.charges == four.charges);
    CHECK(one.gauss_newton == four.gauss_newton);
    CHECK(one.gradient == four.gradient);
    CHECK(one.sum_of_squares == four.sum_of_squares);

    // The rows of every batch are in place.
    const auto charges = charges_of(plan, candidate);
    CHECK(one.charges == charges);
    auto sum_of_squares = 0.0;
    for (std::size_t row = 0; row < charges.size(); ++row) {
        const auto residual = weights[row] * (charges[row] - references[row]);
        sum_of_squares += residual * residual;
    }
    CHECK(std::abs(one.sum_of_squares - sum_of_squares) <= 1e-12 * sum_of_squares);
}

TEST_CASE("Normal equations reject what they cannot calculate", "[calculation][derivatives]") {
    const auto assessment = assess(sqe_molecules(), unusual_sqeqp_parameters());
    const auto& plan = *assessment.default_plan();
    const auto candidate = parameterization::with_parameter_values(
        *plan.parameter_set(), parameterization::parameter_values(*plan.parameter_set()),
        {.id = "candidate"});
    const auto references = std::vector<double>(8, 0.0);
    const auto weights = std::vector<double>(8, 1.0);
    constexpr auto nan = std::numeric_limits<double>::quiet_NaN();

    // Malformed references and weights.
    CHECK_THROWS_AS(
        parameterization::normal_equations(plan, candidate, std::vector<double>(7, 0.0), weights),
        std::invalid_argument);
    auto bad = references;
    bad[2] = nan;
    CHECK_THROWS_AS(parameterization::normal_equations(plan, candidate, bad, weights),
                    std::invalid_argument);
    bad = weights;
    bad[2] = nan;
    CHECK_THROWS_AS(parameterization::normal_equations(plan, candidate, references, bad),
                    std::invalid_argument);
    bad[2] = -1.0;
    CHECK_THROWS_AS(parameterization::normal_equations(plan, candidate, references, bad),
                    std::invalid_argument);

    // The plan's own parameter set ID.
    CHECK_THROWS_AS(
        parameterization::normal_equations(plan, *plan.parameter_set(), references, weights),
        std::invalid_argument);

    // Cutoff execution, and a method without exact derivatives.
    const auto cutoff =
        assess(sqe_molecules(), chargefw::test::make_eem_ho_parameters(0.5), std::nullopt,
               calculation::ExecutionSelection{calculation::ExecutionSelectionKind::cutoff, 8.0});
    const auto& cutoff_plan = *cutoff.default_plan();
    CHECK_THROWS_AS(parameterization::normal_equations(
                        cutoff_plan,
                        parameterization::with_parameter_values(
                            *cutoff_plan.parameter_set(),
                            parameterization::parameter_values(*cutoff_plan.parameter_set()),
                            {.id = "candidate"}),
                        references, weights),
                    std::invalid_argument);
    const auto qeq = assess(sqe_molecules(), chargefw::test::make_qeq_ho_parameters());
    const auto& qeq_plan = *qeq.default_plan();
    CHECK_THROWS_AS(
        parameterization::normal_equations(
            qeq_plan,
            parameterization::with_parameter_values(
                *qeq_plan.parameter_set(),
                parameterization::parameter_values(*qeq_plan.parameter_set()), {.id = "candidate"}),
            references, weights),
        std::invalid_argument);

    // Weights so large that the weighted sums overflow.
    CHECK_THROWS_AS(parameterization::normal_equations(plan, candidate, references,
                                                       std::vector<double>(8, 1e200)),
                    std::runtime_error);
}
