#include "support/test_charge_results.h"

#include <chargefw/adapters/charge_result.h>
#include <chargefw/calculation/calculation.h>
#include <chargefw/charges/atomic_charges.h>
#include <chargefw/charges/charge_collection.h>
#include <chargefw/core/atom.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>

#include <snitch/snitch.hpp>

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace adapters = chargefw::adapters;
namespace calculation = chargefw::calculation;
namespace charges = chargefw::charges;
namespace core = chargefw::core;

using chargefw::test::full_effective;

TEST_CASE("result assembly validates fixed ion provenance", "[adapters][result]") {
    const auto records = std::vector{adapters::ImportedMoleculeRecord{
        .molecule = core::Molecule{std::vector{core::Atom{6}, core::Atom{8}}}}};
    const auto make_result = [&records](calculation::FixedIons fixed_ions) {
        return adapters::make_charge_calculation_result(
            records, {},
            {.status = calculation::ExecutionStatus::cancelled,
             .effective = full_effective("eem", std::nullopt, std::move(fixed_ions))});
    };
    const auto valid =
        calculation::FixedIons{.sources = {{.molecule_index = 0, .atom_index = 1, .charge = 0.5}}};
    CHECK_NOTHROW(make_result(valid));

    const auto source = valid.sources.front();
    const auto invalid_cases = std::vector<calculation::FixedIons>{
        {.sources = {}},
        {.sources = {{.molecule_index = 1, .atom_index = 1, .charge = 0.5}}},
        {.sources = {{.molecule_index = 0, .atom_index = 2, .charge = 0.5}}},
        {.sources = {{.molecule_index = 0,
                      .atom_index = 1,
                      .charge = std::numeric_limits<double>::quiet_NaN()}}},
        {.sources = {source, source}}};
    for (std::size_t index = 0; index < invalid_cases.size(); ++index) {
        CAPTURE(index);
        CHECK_THROWS_AS(make_result(invalid_cases[index]), std::invalid_argument);
    }
}

TEST_CASE("results hold either an execution or a failed request", "[adapters][result]") {
    const auto records = std::vector{adapters::ImportedMoleculeRecord{
        .molecule = core::Molecule{std::vector{core::Atom{1}}},
        .diagnostics = {{.code = "partial_charges_ignored", .message = "ignored"}}}};

    const auto cancelled = adapters::make_charge_calculation_result(
        records, {}, {.status = calculation::ExecutionStatus::cancelled});
    CHECK(cancelled.status() == adapters::ResultStatus::cancelled);
    REQUIRE(cancelled.execution() != nullptr);
    CHECK_FALSE(cancelled.request_failure().has_value());
    CHECK(cancelled.charges() == nullptr);

    const auto failed = adapters::make_failed_request_result(records, {}, "bad threshold");
    CHECK(failed.status() == adapters::ResultStatus::invalid_request);
    CHECK(failed.execution() == nullptr);
    CHECK(failed.request_failure() == "bad threshold");
    CHECK(failed.charges() == nullptr);
    CHECK(adapters::to_string(failed.status()) == "invalid_input_or_request");
    // The record keeps its import warnings and reports the failure after them.
    const auto diagnostics = adapters::charge_record_diagnostics(failed, 0);
    REQUIRE(diagnostics.size() == 2);
    CHECK(diagnostics[0].code == "partial_charges_ignored");
    CHECK(diagnostics[1].code == "invalid_input_or_request");
    CHECK(diagnostics[1].message == "bad threshold");
}

TEST_CASE("record diagnostics carry the input warnings of their own molecule",
          "[adapters][result]") {
    const auto records = std::vector(
        3,
        adapters::ImportedMoleculeRecord{.molecule = core::Molecule{std::vector{core::Atom{6}}}});
    auto execution =
        calculation::ExecutionResult{.status = calculation::ExecutionStatus::cancelled};
    for (const auto molecule_index : {1U, 2U}) {
        execution.input_warnings.push_back({.kind = calculation::InputIssueKind::hydrogens_missing,
                                            .message = "no hydrogens",
                                            .molecule_index = molecule_index});
    }
    const auto result = adapters::make_charge_calculation_result(records, {}, std::move(execution));

    for (std::size_t molecule_index = 0; molecule_index < records.size(); ++molecule_index) {
        CAPTURE(molecule_index);
        auto warnings = std::vector<adapters::ResultDiagnostic>{};
        std::ranges::copy_if(adapters::charge_record_diagnostics(result, molecule_index),
                             std::back_inserter(warnings), [](const auto& diagnostic) {
                                 return diagnostic.code == "hydrogens_missing";
                             });
        REQUIRE(warnings.size() == (molecule_index == 0 ? 0U : 1U));
        if (!warnings.empty()) {
            CHECK(warnings[0].severity == adapters::DiagnosticSeverity::warning);
            CHECK(warnings[0].molecule_index == molecule_index);
        }
    }

    // Warnings out of molecule order, or for a molecule outside the input, are refused.
    for (const auto& molecule_indices :
         {std::vector<std::size_t>{2, 1}, std::vector<std::size_t>{3}}) {
        auto invalid =
            calculation::ExecutionResult{.status = calculation::ExecutionStatus::cancelled};
        for (const auto molecule_index : molecule_indices) {
            invalid.input_warnings.push_back(
                {.kind = calculation::InputIssueKind::hydrogens_missing,
                 .message = "no hydrogens",
                 .molecule_index = molecule_index});
        }
        CHECK_THROWS_AS(adapters::make_charge_calculation_result(records, {}, std::move(invalid)),
                        std::invalid_argument);
    }
}

TEST_CASE("result assembly validates assignment dimensions targets and scope",
          "[adapters][result]") {
    const auto records = std::vector{adapters::ImportedMoleculeRecord{
        .molecule = core::Molecule{std::vector{core::Atom{1}, core::Atom{1}},
                                   {},
                                   std::vector{core::Conformer{{core::Position{0.0, 0.0, 0.0},
                                                                core::Position{1.0, 0.0, 0.0}}}},
                                   "hydrogen"},
        .identity = {.source = "hydrogen.json", .record_id = "hydrogen"}}};
    const auto pair_charges = charges::AtomicCharges{{0.0, 0.0}};
    const auto formal = [](std::vector<charges::ChargeAssignment> assignments) {
        return charges::ChargeSet{"formal", std::move(assignments)};
    };

    // Successful results without effective provenance receive full-execution provenance, so each
    // case isolates its own invalid charge data or status.
    auto invalid_cases = std::vector<calculation::ExecutionResult>{};
    invalid_cases.push_back({});
    invalid_cases.push_back({.charges = formal({{.target = {.molecule_index = 0},
                                                 .charges = charges::AtomicCharges{{0.0}}}})});
    invalid_cases.push_back(
        {.charges = formal({{.target = {.molecule_index = 1}, .charges = pair_charges}})});
    invalid_cases.push_back(
        {.charges = formal(
             {{.target = {.molecule_index = 0, .conformer_index = 1}, .charges = pair_charges}})});
    invalid_cases.push_back(
        {.charges = formal(
             {{.target = {.molecule_index = 0}, .charges = pair_charges},
              {.target = {.molecule_index = 0, .conformer_index = 0}, .charges = pair_charges}})});
    invalid_cases.push_back(
        {.status = calculation::ExecutionStatus::numerical_failure,
         .charges = formal({{.target = {.molecule_index = 0}, .charges = pair_charges}})});

    for (std::size_t index = 0; index < invalid_cases.size(); ++index) {
        CAPTURE(index);
        auto result = std::move(invalid_cases[index]);
        if (result.status == calculation::ExecutionStatus::success &&
            !result.effective.has_value()) {
            result.effective = full_effective("formal");
        }
        CHECK_THROWS_AS(adapters::make_charge_calculation_result(records, {}, std::move(result)),
                        std::invalid_argument);
    }

    auto invalid_mapping = records;
    invalid_mapping[0].import_metadata = adapters::MoleculeImportMetadata{
        .format = adapters::MolecularSourceFormat::molecule_json,
        .atoms = {{.position = 0}},
        .conformers = {{.position = 0, .sites = {{.position = 0}}}}};
    CHECK_THROWS_AS(
        adapters::make_charge_calculation_result(
            std::move(invalid_mapping), {},
            {.charges = formal({{.target = {.molecule_index = 0}, .charges = pair_charges}}),
             .effective = full_effective("formal")}),
        std::invalid_argument);
}

TEST_CASE("result assembly requires canonical assignment order", "[adapters][result]") {
    const auto records =
        std::vector{adapters::ImportedMoleculeRecord{.molecule = core::Molecule{{core::Atom{1}}}},
                    adapters::ImportedMoleculeRecord{.molecule = core::Molecule{{core::Atom{1}}}}};
    const auto make_result = [&records](const std::size_t first, const std::size_t second) {
        return adapters::make_charge_calculation_result(
            records, {},
            {.charges = charges::ChargeSet{"formal",
                                           {{.target = {.molecule_index = first},
                                             .charges = charges::AtomicCharges{{0.0}}},
                                            {.target = {.molecule_index = second},
                                             .charges = charges::AtomicCharges{{0.0}}}}},
             .effective = full_effective("formal")});
    };

    CHECK_NOTHROW(make_result(0, 1));
    CHECK_THROWS_AS(make_result(1, 0), std::invalid_argument);
}
