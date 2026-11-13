#pragma once

#include "support/test_molecules.h"

#include <chargefw/calculation/calculation.h>
#include <chargefw/calculation/observer.h>
#include <chargefw/charges/charge_set.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/parameters/models/parameter_set.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <snitch/snitch.hpp>

namespace chargefw::test {

/// Assesses the request and calculates its default plan.
[[nodiscard]] inline auto calculate_application(
    calculation::AssessmentRequest request,
    const calculation::CalculationObserver& observer = calculation::default_calculation_observer(),
    const std::size_t max_threads = 1) -> calculation::ExecutionResult {
    const auto assessment = calculation::assess(std::move(request));
    return calculation::calculate(assessment, max_threads, observer);
}

[[nodiscard]] inline auto
calculate_method(core::Molecule molecule, std::string_view method_id,
                 std::vector<parameters::ParameterSet> parameter_sets = {},
                 const methods::MethodOptions* method_options = nullptr) -> charges::ChargeSet {
    auto request = calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{std::move(molecule)}, "test"},
        .parameter_sets = std::move(parameter_sets),
        .method_id = std::string{method_id}};

    if (method_options != nullptr) {
        request.method_options.emplace(method_id, *method_options);
    }

    auto assessment = calculation::assess(std::move(request));
    REQUIRE(assessment.default_plan() != nullptr);

    auto result = calculation::calculate(assessment);
    REQUIRE(result.calculated());
    return std::move(*result.charges);
}

[[nodiscard]] inline auto
calculate_single_method(core::Molecule molecule, std::string_view method_id,
                        std::vector<parameters::ParameterSet> parameter_sets = {},
                        const methods::MethodOptions* method_options = nullptr)
    -> charges::ChargeSet {
    const auto charge_set =
        calculate_method(std::move(molecule), method_id, std::move(parameter_sets), method_options);
    REQUIRE(charge_set.size() == 1);
    return charge_set;
}

inline auto assert_conformer_dependent(const charges::ChargeSet& charge_set,
                                       const std::size_t conformer_count) -> void {
    REQUIRE(charge_set.size() == conformer_count);

    for (std::size_t conformer_index = 0; conformer_index < conformer_count; ++conformer_index) {
        const auto& assignment = charge_set.assignment(conformer_index);
        CHECK(assignment.target.molecule_index == 0);
        CHECK(assignment.target.conformer_index == conformer_index);
    }
}

/// Asserts the chemically expected charge shape on neutral water: the oxygen carries the
/// negative charge, both hydrogens carry equal positive charge, and the molecule is neutral.
inline auto assert_neutral_water_charges(const charges::AtomicCharges& charges,
                                         const double tolerance) -> void {
    REQUIRE(charges.size() == 3);
    CHECK(charges[0] < 0.0);
    CHECK(charges[1] > 0.0);
    CHECK(charges[2] > 0.0);
    CHECK(std::abs(charges[1] - charges[2]) < tolerance);
    CHECK(std::abs(charges.total() - 0.0) < tolerance);
}

inline auto assert_same_charges(const charges::AtomicCharges& actual,
                                const charges::AtomicCharges& expected, const double tolerance)
    -> void {
    REQUIRE(actual.size() == expected.size());

    for (std::size_t atom_index = 0; atom_index < actual.size(); ++atom_index) {
        CHECK(std::abs(actual[atom_index] - expected[atom_index]) < tolerance);
    }
}

/// Method-agnostic invariant: charges must not depend on atom numbering or bond endpoint
/// order. Verifies a relabeled water (same physical molecule, permuted indices) and a
/// bond-flipped water against the plain-water reference.
inline auto
assert_water_charges_labeling_invariant(std::string_view method_id,
                                        const charges::AtomicCharges& reference_charges,
                                        std::vector<parameters::ParameterSet> parameter_sets = {},
                                        const methods::MethodOptions* method_options = nullptr,
                                        const double tolerance = 1.0e-9) -> void {
    static constexpr std::array new_to_old{std::size_t{2}, std::size_t{0}, std::size_t{1}};

    const auto relabeled = calculate_single_method(relabel_atoms(make_water(), new_to_old),
                                                   method_id, parameter_sets, method_options);
    const auto& relabeled_charges = relabeled.assignment(0).charges;

    REQUIRE(relabeled_charges.size() == reference_charges.size());
    for (std::size_t atom_index = 0; atom_index < new_to_old.size(); ++atom_index) {
        CHECK(std::abs(relabeled_charges[atom_index] -
                       reference_charges[new_to_old.at(atom_index)]) < tolerance);
    }

    const auto flipped = calculate_single_method(flip_bond_directions(make_water()), method_id,
                                                 std::move(parameter_sets), method_options);
    assert_same_charges(flipped.assignment(0).charges, reference_charges, tolerance);
}

} // namespace chargefw::test
