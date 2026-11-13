#include "support/test_charge_results.h"

#include <chargefw/adapters/charge_result.h>
#include <chargefw/calculation/calculation.h>
#include <chargefw/charges/atomic_charges.h>
#include <chargefw/charges/charge_set.h>
#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/methods/method_prerequisites.h>

#include <snitch/snitch.hpp>

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace adapters = chargefw::adapters;
namespace calculation = chargefw::calculation;
namespace charges = chargefw::charges;
namespace core = chargefw::core;
namespace methods = chargefw::methods;

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
        .issues = {
            {.kind = adapters::ImportIssueKind::partial_charges_ignored, .message = "ignored"}}}};

    const auto cancelled = adapters::make_charge_calculation_result(
        records, {}, {.status = calculation::ExecutionStatus::cancelled});
    CHECK(cancelled.status() == adapters::ResultStatus::cancelled);
    REQUIRE(cancelled.execution() != nullptr);
    CHECK_FALSE(cancelled.request_failure().has_value());
    CHECK(cancelled.charges() == nullptr);

    const auto failed = adapters::make_failed_request_result(records, {}, "bad threshold");
    CHECK(failed.status() == adapters::ResultStatus::invalid_input_or_request);
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
    execution.input_warnings.push_back({.kind = calculation::InputIssueKind::hydrogens_missing,
                                        .message = "no hydrogens",
                                        .molecule_index = 1});
    execution.input_warnings.push_back({.kind = calculation::InputIssueKind::atoms_too_close,
                                        .message = "atoms too close",
                                        .molecule_index = 2,
                                        .location = {.conformer_index = 0, .atom_index = 0}});
    const auto result = adapters::make_charge_calculation_result(records, {}, std::move(execution));

    for (std::size_t molecule_index = 0; molecule_index < records.size(); ++molecule_index) {
        CAPTURE(molecule_index);
        auto warnings = std::vector<adapters::ResultDiagnostic>{};
        std::ranges::copy_if(adapters::charge_record_diagnostics(result, molecule_index),
                             std::back_inserter(warnings), [](const auto& diagnostic) {
                                 return diagnostic.code == "hydrogens_missing" ||
                                        diagnostic.code == "atoms_too_close";
                             });
        REQUIRE(warnings.size() == (molecule_index == 0 ? 0U : 1U));
        if (!warnings.empty()) {
            CHECK(warnings[0].severity == adapters::DiagnosticSeverity::warning);
            CHECK(warnings[0].molecule_index == molecule_index);
        }
    }
    const auto close = adapters::charge_record_diagnostics(result, 2).front();
    CHECK(close.location.conformer_index == 0);
    CHECK(close.location.atom_index == 0);

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

TEST_CASE("record diagnostics carry the charge warnings with all their atoms",
          "[adapters][result]") {
    const auto records =
        std::vector(2, adapters::ImportedMoleculeRecord{
                           .molecule = core::Molecule{std::vector{core::Atom{6}, core::Atom{8}}}});
    const auto make_execution = [](const calculation::ExecutionStatus status) {
        auto execution = calculation::ExecutionResult{.status = status};
        if (status == calculation::ExecutionStatus::success) {
            execution.charges = charges::ChargeSet{
                "dummy",
                {{.target = {.molecule_index = 0}, .charges = charges::AtomicCharges{{0.0, 0.0}}},
                 {.target = {.molecule_index = 1}, .charges = charges::AtomicCharges{{0.0, 0.0}}}}};
            execution.effective = full_effective("dummy");
        }
        execution.charge_warnings.push_back(
            {.kind = calculation::ChargeIssueKind::charge_out_of_range,
             .message = "charges out of range",
             .molecule_index = 1,
             .location = {.atom_index = 1, .atom_indices = {0, 1}}});
        // The same text for other atoms is another finding.
        execution.charge_warnings.push_back(
            {.kind = calculation::ChargeIssueKind::charge_out_of_range,
             .message = "charges out of range",
             .molecule_index = 1,
             .location = {.atom_index = 0, .atom_indices = {0}}});
        return execution;
    };
    const auto result = adapters::make_charge_calculation_result(
        records, {}, make_execution(calculation::ExecutionStatus::success));

    CHECK(adapters::charge_record_diagnostics(result, 0).empty());
    const auto diagnostics = adapters::charge_record_diagnostics(result, 1);
    REQUIRE(diagnostics.size() == 2);
    CHECK(diagnostics[0].code == "charge_out_of_range");
    CHECK(diagnostics[0].severity == adapters::DiagnosticSeverity::warning);
    CHECK(diagnostics[0].location.atom_index == 1);
    CHECK(diagnostics[0].location.atom_indices == std::vector<std::size_t>{0, 1});
    CHECK(diagnostics[1].location.atom_indices == std::vector<std::size_t>{0});

    // Charge warnings without charges are refused.
    CHECK_THROWS_AS(adapters::make_charge_calculation_result(
                        records, {}, make_execution(calculation::ExecutionStatus::cancelled)),
                    std::invalid_argument);
}

TEST_CASE("record diagnostics name the residues of the atoms they point to", "[adapters][result]") {
    const auto labels = [](std::string residue, std::string chain, std::string sequence) {
        auto value = adapters::SourceStructuralLabels{};
        value.author = {.residue = std::move(residue),
                        .chain = std::move(chain),
                        .sequence = std::move(sequence)};
        value.insertion_code = "?";
        return std::make_shared<const adapters::SourceStructuralLabels>(std::move(value));
    };
    const auto molecule = core::Molecule{std::vector{core::Atom{6}, core::Atom{8}, core::Atom{8}},
                                         std::vector{core::Bond{0, 1}, core::Bond{1, 2}}};
    auto structure = adapters::ImportedMoleculeRecord{.molecule = molecule};
    structure.import_metadata = adapters::MoleculeImportMetadata{
        .format = adapters::MolecularSourceFormat::mmcif,
        .atoms = {{.position = 0, .structural_labels = labels("ALA", "A", "5")},
                  {.position = 1, .structural_labels = labels("HOH", "A", "1203")},
                  {.position = 2, .structural_labels = labels("HOH", "A", "1203")}}};
    // A molecule built in code has no residues, so its diagnostics keep the plain message.
    const auto records =
        std::vector{structure, adapters::ImportedMoleculeRecord{.molecule = molecule}};

    auto execution =
        calculation::ExecutionResult{.status = calculation::ExecutionStatus::no_executable_plan};
    for (const auto molecule_index : {0U, 1U}) {
        execution.rejections.push_back(
            {.method_id = "sqe",
             .issues = {methods::PrerequisiteIssue{
                            .kind = methods::PrerequisiteIssueKind::parameter_classification_failed,
                            .message = "no atom parameters",
                            .molecule_index = molecule_index,
                            .location = {.atom_index = 1, .atom_indices = {1, 2}}},
                        methods::PrerequisiteIssue{
                            .kind = methods::PrerequisiteIssueKind::parameter_classification_failed,
                            .message = "no bond parameters",
                            .molecule_index = molecule_index,
                            .location = {.bond_index = 0, .bond_indices = {0}}},
                        methods::PrerequisiteIssue{
                            .kind = methods::PrerequisiteIssueKind::invalid_geometry,
                            .message = "identical coordinates",
                            .molecule_index = molecule_index,
                            .location = {.atom_index = 2, .other_atom_index = 1}}}});
    }
    const auto result = adapters::make_charge_calculation_result(records, {}, std::move(execution));

    const auto messages = [&result](const std::size_t molecule_index) {
        auto values = std::vector<std::string>{};
        for (const auto& diagnostic : adapters::charge_record_diagnostics(result, molecule_index)) {
            if (diagnostic.severity == adapters::DiagnosticSeverity::error &&
                diagnostic.code != "no_executable_plan") {
                values.push_back(diagnostic.message);
            }
        }
        return values;
    };
    // Residues are only added: each message starts with the one for the molecule without them.
    const auto plain = messages(1);
    const auto enriched = messages(0);
    REQUIRE(plain.size() == 3);
    REQUIRE(enriched.size() == 3);
    auto added = std::vector<std::string>{};
    for (std::size_t index = 0; index < plain.size(); ++index) {
        CAPTURE(plain[index], enriched[index]);
        CHECK_FALSE(plain[index].contains("HOH"));
        REQUIRE(enriched[index].starts_with(plain[index]));
        added.push_back(enriched[index].substr(plain[index].size()));
        // The placeholder insertion code is left out.
        CHECK_FALSE(added.back().contains('?'));
    }
    // The atom, both atoms of the bond, and both atoms of the pair, with their residues.
    CHECK(added[0].contains("atom 2"));
    CHECK(added[0].contains("HOH A 1203"));
    CHECK_FALSE(added[0].contains("ALA"));
    CHECK(added[1].contains("ALA A 5"));
    CHECK(added[1].contains("HOH A 1203"));
    CHECK(added[2].contains("atoms 2 and 3"));
    CHECK(added[2].contains("HOH A 1203"));
}

TEST_CASE("plan warnings go to their molecule's record or to the document", "[adapters][result]") {
    const auto records = std::vector(
        2,
        adapters::ImportedMoleculeRecord{.molecule = core::Molecule{std::vector{core::Atom{6}}}});
    const auto make_result = [&records](
                                 std::vector<chargefw::calculation::ExecutionIssue> warnings) {
        auto effective = full_effective("eem");
        effective.warnings = std::move(warnings);
        return adapters::make_charge_calculation_result(
            records, {},
            {.status = calculation::ExecutionStatus::cancelled, .effective = std::move(effective)});
    };
    const auto threshold = [](std::optional<std::size_t> molecule_index) {
        return chargefw::calculation::ExecutionIssue{
            .kind = chargefw::calculation::ExecutionIssueKind::resource_threshold_exceeded,
            .message = "above the limit",
            .molecule_index = molecule_index};
    };
    const auto codes = [](const std::vector<adapters::ResultDiagnostic>& diagnostics) {
        auto result = std::vector<std::string>{};
        for (const auto& diagnostic : diagnostics) {
            result.push_back(diagnostic.code);
        }
        return result;
    };

    const auto result = make_result({threshold(std::nullopt), threshold(1)});
    CHECK(codes(adapters::charge_result_diagnostics(result)) ==
          std::vector<std::string>{"resource_threshold_exceeded", "calculation_cancelled"});
    CHECK(codes(adapters::charge_record_diagnostics(result, 0)) ==
          std::vector<std::string>{"calculation_cancelled"});
    CHECK(codes(adapters::charge_record_diagnostics(result, 1)) ==
          std::vector<std::string>{"resource_threshold_exceeded", "calculation_cancelled"});

    CHECK_THROWS_AS(make_result({threshold(1), threshold(0)}), std::invalid_argument);
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
