#include <chargefw/adapters/charge_result.h>

#include "adapters/record_assignments.h"

#include <chargefw/calculation/calculation.h>
#include <chargefw/methods/method_prerequisites.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iterator>
#include <memory>
#include <optional>
#include <ranges>
#include <set>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

namespace chargefw::adapters {
namespace {

[[nodiscard]] auto prerequisite_code(const methods::PrerequisiteIssueKind kind)
    -> std::string_view {
    switch (kind) {
    case methods::PrerequisiteIssueKind::invalid_options:
        return "invalid_method_options";
    case methods::PrerequisiteIssueKind::missing_feature:
        return "missing_feature";
    case methods::PrerequisiteIssueKind::invalid_geometry:
        return "invalid_geometry";
    case methods::PrerequisiteIssueKind::unsupported_molecule:
        return "unsupported_molecule";
    case methods::PrerequisiteIssueKind::missing_parameters:
        return "missing_parameters";
    case methods::PrerequisiteIssueKind::parameter_classification_failed:
        return "parameter_classification_failed";
    case methods::PrerequisiteIssueKind::unsupported_fixed_ions:
        return "unsupported_fixed_ions";
    }
    throw std::logic_error{"unknown prerequisite issue kind"};
}

// Code and message of the diagnostics already appended; a molecule can have one per atom.
using SeenDiagnostics = std::set<std::pair<std::string, std::string>>;

auto append_unique(std::vector<ResultDiagnostic>& diagnostics, SeenDiagnostics& seen,
                   ResultDiagnostic diagnostic) -> void {
    if (seen.emplace(diagnostic.code, diagnostic.message).second) {
        diagnostics.push_back(std::move(diagnostic));
    }
}

auto validate_import_metadata(const std::span<const ImportedMoleculeRecord> records) -> void {
    for (const auto& record : records) {
        if (record.caller_atom_ids.has_value() &&
            record.caller_atom_ids->size() != record.molecule.atom_count()) {
            throw std::invalid_argument{"caller atom ID count does not match molecule atom count"};
        }
        if (record.caller_atom_ids.has_value() && record.import_metadata.has_value()) {
            throw std::invalid_argument{"caller atom IDs cannot accompany import metadata"};
        }
        if (!record.import_metadata.has_value()) {
            continue;
        }
        const auto& metadata = *record.import_metadata;
        if (metadata.atoms.size() != record.molecule.atom_count()) {
            throw std::invalid_argument{
                "import atom mapping size does not match molecule atom count"};
        }
        if (metadata.conformers.size() != record.molecule.conformer_count()) {
            throw std::invalid_argument{
                "import conformer mapping size does not match molecule conformer count"};
        }
        for (const auto& conformer : metadata.conformers) {
            if (conformer.sites.size() != record.molecule.atom_count()) {
                throw std::invalid_argument{
                    "import conformer site mapping size does not match molecule atom count"};
            }
        }
    }
}

auto validate_fixed_ions(const std::span<const ImportedMoleculeRecord> records,
                         const calculation::EffectiveCalculation& effective) -> void {
    if (!effective.fixed_ions.has_value()) {
        return;
    }

    const auto& fixed_ions = *effective.fixed_ions;
    if (fixed_ions.sources.empty()) {
        throw std::invalid_argument{"fixed ion provenance requires sources"};
    }
    auto selectors = std::set<std::pair<std::size_t, std::size_t>>{};
    for (const auto& source : fixed_ions.sources) {
        if (source.molecule_index >= records.size()) {
            throw std::invalid_argument{"fixed ion source molecule index is outside "
                                        "the input"};
        }
        if (source.atom_index >= records[source.molecule_index].molecule.atom_count()) {
            throw std::invalid_argument{"fixed ion source atom index is outside its "
                                        "molecule"};
        }
        if (!std::isfinite(source.charge)) {
            throw std::invalid_argument{"fixed ion source charge must be finite"};
        }
        if (!selectors.emplace(source.molecule_index, source.atom_index).second) {
            throw std::invalid_argument{"fixed ion source selectors must be unique"};
        }
    }
}

// The warnings of the plan that ran. Record diagnostics find a molecule's warnings by binary
// search.
[[nodiscard]] auto plan_warnings(const calculation::ExecutionResult& execution)
    -> std::span<const methods::ExecutionIssue> {
    if (!execution.effective.has_value()) {
        return {};
    }
    return execution.effective->execution_issues;
}

auto validate_plan_warnings(const std::span<const ImportedMoleculeRecord> records,
                            const calculation::ExecutionResult& execution) -> void {
    const auto warnings = plan_warnings(execution);
    if (!std::ranges::is_sorted(warnings, {}, &methods::ExecutionIssue::molecule_index)) {
        throw std::invalid_argument{"plan warnings are not in molecule order"};
    }
    if (!warnings.empty() && warnings.back().molecule_index.has_value() &&
        *warnings.back().molecule_index >= records.size()) {
        throw std::invalid_argument{"plan warning molecule index is outside the input"};
    }
}

[[nodiscard]] auto plan_warning_diagnostic(const methods::ExecutionIssue& issue)
    -> ResultDiagnostic {
    auto diagnostic = ResultDiagnostic{};
    diagnostic.severity = DiagnosticSeverity::warning;
    diagnostic.code = methods::to_string(issue.kind);
    diagnostic.message = issue.message;
    diagnostic.molecule_index = issue.molecule_index;
    return diagnostic;
}

// Record diagnostics find a molecule's input warnings by binary search.
auto validate_input_warnings(const std::span<const ImportedMoleculeRecord> records,
                             const calculation::ExecutionResult& execution) -> void {
    const auto& warnings = execution.input_warnings;
    if (!std::ranges::is_sorted(warnings, {}, &calculation::InputIssue::molecule_index)) {
        throw std::invalid_argument{"input warnings are not in molecule order"};
    }
    if (!warnings.empty() && warnings.back().molecule_index >= records.size()) {
        throw std::invalid_argument{"input warning molecule index is outside the input"};
    }
}

auto validate_assignments(const std::span<const ImportedMoleculeRecord> records,
                          const calculation::ExecutionResult& result) -> void {
    if (result.effective.has_value()) {
        validate_fixed_ions(records, *result.effective);
    }
    if (result.status != calculation::ExecutionStatus::success) {
        if (!result.charges.has_value()) {
            return;
        }
        throw std::invalid_argument{"unsuccessful calculation result contains charge assignments"};
    }
    if (!result.charges.has_value()) {
        throw std::invalid_argument{"successful calculation result has no charge assignments"};
    }
    if (!result.effective.has_value()) {
        throw std::invalid_argument{"successful calculation result has no effective provenance"};
    }
    if (result.charges->method_id() != result.effective->method_id) {
        throw std::invalid_argument{"charge assignment method does not match effective provenance"};
    }
    const auto charge_parameter_set = result.charges->parameter_set_id();
    if (charge_parameter_set.has_value() != result.effective->parameter_set_id.has_value() ||
        (charge_parameter_set.has_value() &&
         *charge_parameter_set != *result.effective->parameter_set_id)) {
        throw std::invalid_argument{
            "charge assignment parameter set does not match effective provenance"};
    }

    const auto assignments = result.charges->assignments();
    if (assignments.empty()) {
        if (records.empty()) {
            return;
        }
        throw std::invalid_argument{"calculation result does not cover every input molecule"};
    }
    const auto conformer_scope = assignments.front().target.conformer_index.has_value();
    const auto validate_assignment = [&](const charges::ChargeAssignment& assignment) {
        const auto molecule_index = assignment.target.molecule_index;
        if (molecule_index >= records.size()) {
            throw std::invalid_argument{"charge assignment molecule index is outside the input"};
        }
        const auto& molecule = records[molecule_index].molecule;
        if (assignment.charges.size() != molecule.atom_count()) {
            throw std::invalid_argument{
                "charge assignment size does not match molecule atom count"};
        }
        if (assignment.target.conformer_index.has_value() != conformer_scope) {
            throw std::invalid_argument{
                "calculation result mixes molecule and conformer assignment scopes"};
        }
        if (assignment.target.conformer_index.has_value() &&
            *assignment.target.conformer_index >= molecule.conformer_count()) {
            throw std::invalid_argument{
                "charge assignment conformer index is outside the molecule"};
        }
    };

    auto assignment_index = std::size_t{0};
    const auto consume = [&](const std::size_t molecule_index,
                             const std::optional<std::size_t> conformer_index) {
        if (assignment_index == assignments.size()) {
            throw std::invalid_argument{"calculation result does not cover every input molecule"};
        }
        const auto& assignment = assignments[assignment_index];
        validate_assignment(assignment);
        if (assignment.target.molecule_index != molecule_index ||
            assignment.target.conformer_index != conformer_index) {
            throw std::invalid_argument{"charge assignments are not in canonical input order"};
        }
        ++assignment_index;
    };
    for (std::size_t molecule_index = 0; molecule_index < records.size(); ++molecule_index) {
        const auto& molecule = records[molecule_index].molecule;
        if (conformer_scope) {
            if (molecule.conformer_count() == 0) {
                throw std::invalid_argument{
                    "calculation result does not cover every input molecule"};
            }
            for (std::size_t conformer_index = 0; conformer_index < molecule.conformer_count();
                 ++conformer_index) {
                consume(molecule_index, conformer_index);
            }
        } else {
            consume(molecule_index, std::nullopt);
        }
    }
    if (assignment_index != assignments.size()) {
        validate_assignment(assignments[assignment_index]);
        throw std::invalid_argument{"charge assignments are not in canonical input order"};
    }
}

[[nodiscard]] auto status_diagnostic(const ChargeCalculationResult& result)
    -> std::optional<ResultDiagnostic> {
    const auto make = [](const DiagnosticSeverity severity, std::string code, std::string message) {
        auto diagnostic = ResultDiagnostic{};
        diagnostic.severity = severity;
        diagnostic.code = std::move(code);
        diagnostic.message = std::move(message);
        return diagnostic;
    };
    const auto* execution = result.execution();
    const auto failure_message = [execution](std::string fallback) {
        return execution->failure_message.value_or(std::move(fallback));
    };
    switch (result.status()) {
    case ResultStatus::success:
        return std::nullopt;
    case ResultStatus::invalid_request:
        return make(DiagnosticSeverity::error,
                    std::string{to_string(ResultStatus::invalid_request)},
                    std::string{*result.request_failure()});
    case ResultStatus::no_executable_plan:
        return make(DiagnosticSeverity::error, "no_executable_plan",
                    "No executable method and parameter-set plan was found.");
    case ResultStatus::numerical_failure:
        return make(DiagnosticSeverity::error, "numerical_failure",
                    failure_message("Calculation failed numerically."));
    case ResultStatus::cancelled:
        return make(DiagnosticSeverity::info, "calculation_cancelled",
                    "Calculation was cancelled before completion.");
    }
    throw std::logic_error{"unknown result status"};
}

[[nodiscard]] auto issue_molecule_index(const calculation::RejectionIssue& issue_value)
    -> std::optional<std::size_t> {
    if (const auto* issue = std::get_if<methods::PrerequisiteIssue>(&issue_value)) {
        return issue->molecule_index;
    }
    return std::get<methods::ExecutionIssue>(issue_value).molecule_index;
}

} // namespace

auto make_charge_calculation_result(std::vector<ImportedMoleculeRecord> inputs,
                                    RequestedCalculationProvenance requested,
                                    calculation::ExecutionResult execution)
    -> ChargeCalculationResult {
    validate_import_metadata(inputs);
    validate_assignments(inputs, execution);
    validate_input_warnings(inputs, execution);
    validate_plan_warnings(inputs, execution);
    auto result = ChargeCalculationResult{};
    result.inputs_ = std::move(inputs);
    result.requested_ = std::move(requested);
    const auto& stored =
        result.outcome_.emplace<calculation::ExecutionResult>(std::move(execution));
    if (stored.status == calculation::ExecutionStatus::no_executable_plan) {
        result.molecule_issues_.resize(result.inputs_.size());
        const auto& rejections = stored.rejections;
        for (std::size_t rejection_index = 0; rejection_index < rejections.size();
             ++rejection_index) {
            const auto& issues = rejections[rejection_index].issues;
            for (std::size_t issue_index = 0; issue_index < issues.size(); ++issue_index) {
                const auto position = std::pair{rejection_index, issue_index};
                const auto molecule_index = issue_molecule_index(issues[issue_index]);
                if (!molecule_index.has_value()) {
                    result.shared_issues_.push_back(position);
                } else if (*molecule_index < result.molecule_issues_.size()) {
                    result.molecule_issues_[*molecule_index].push_back(position);
                }
            }
        }
    }
    return result;
}

auto ChargeCalculationResult::inputs() const noexcept -> std::span<const ImportedMoleculeRecord> {
    return inputs_;
}

auto ChargeCalculationResult::requested() const noexcept -> const RequestedCalculationProvenance& {
    return requested_;
}

auto make_failed_request_result(std::vector<ImportedMoleculeRecord> inputs,
                                RequestedCalculationProvenance requested, std::string message)
    -> ChargeCalculationResult {
    validate_import_metadata(inputs);
    auto result = ChargeCalculationResult{};
    result.inputs_ = std::move(inputs);
    result.requested_ = std::move(requested);
    result.outcome_.emplace<ChargeCalculationResult::RequestFailure>(std::move(message));
    return result;
}

auto to_string(const ResultStatus value) -> std::string_view {
    switch (value) {
    case ResultStatus::success:
        return "success";
    case ResultStatus::invalid_request:
        return "invalid_input_or_request";
    case ResultStatus::no_executable_plan:
        return "no_executable_plan";
    case ResultStatus::numerical_failure:
        return "numerical_failure";
    case ResultStatus::cancelled:
        return "cancelled";
    }
    throw std::invalid_argument{"unknown result status"};
}

auto ChargeCalculationResult::status() const noexcept -> ResultStatus {
    const auto* execution = this->execution();
    if (execution == nullptr) {
        return ResultStatus::invalid_request;
    }
    switch (execution->status) {
    case calculation::ExecutionStatus::success:
        return ResultStatus::success;
    case calculation::ExecutionStatus::no_executable_plan:
        return ResultStatus::no_executable_plan;
    case calculation::ExecutionStatus::numerical_failure:
        return ResultStatus::numerical_failure;
    case calculation::ExecutionStatus::cancelled:
        return ResultStatus::cancelled;
    }
    std::abort();
}

auto ChargeCalculationResult::execution() const noexcept -> const calculation::ExecutionResult* {
    return std::get_if<calculation::ExecutionResult>(&outcome_);
}

auto ChargeCalculationResult::request_failure() const noexcept -> std::optional<std::string_view> {
    if (const auto* failure = std::get_if<RequestFailure>(&outcome_)) {
        return failure->message;
    }
    return std::nullopt;
}

auto ChargeCalculationResult::charges() const noexcept -> const charges::ChargeSet* {
    const auto* execution = this->execution();
    return execution != nullptr && execution->calculated() ? std::addressof(*execution->charges)
                                                           : nullptr;
}

auto charge_result_diagnostics(const ChargeCalculationResult& result)
    -> std::vector<ResultDiagnostic> {
    auto diagnostics = std::vector<ResultDiagnostic>{};
    if (const auto* execution = result.execution()) {
        // Warnings without a molecule come first in molecule order.
        for (const auto& issue :
             std::ranges::equal_range(plan_warnings(*execution), std::optional<std::size_t>{}, {},
                                      &methods::ExecutionIssue::molecule_index)) {
            diagnostics.push_back(plan_warning_diagnostic(issue));
        }
    }
    if (const auto diagnostic = status_diagnostic(result); diagnostic.has_value()) {
        diagnostics.push_back(*diagnostic);
    }
    return diagnostics;
}

auto charge_record_diagnostics(const ChargeCalculationResult& result,
                               const std::size_t molecule_index) -> std::vector<ResultDiagnostic> {
    const auto records = result.inputs();
    if (molecule_index >= records.size()) {
        throw std::out_of_range{"result diagnostic molecule index is outside the input"};
    }
    const auto* execution = result.execution();
    auto diagnostics = std::vector<ResultDiagnostic>{};
    auto seen = SeenDiagnostics{};
    for (const auto& import_diagnostic : records[molecule_index].diagnostics) {
        auto imported = ResultDiagnostic{};
        imported.severity = DiagnosticSeverity::warning;
        imported.code = import_diagnostic.code;
        imported.message = import_diagnostic.message;
        imported.molecule_index = molecule_index;
        imported.line = import_diagnostic.line;
        append_unique(diagnostics, seen, std::move(imported));
    }
    // A request that failed has no execution result, so only its status follows.
    const auto input_warnings =
        execution != nullptr ? std::span<const calculation::InputIssue>{execution->input_warnings}
                             : std::span<const calculation::InputIssue>{};
    // Input warnings are in molecule order.
    for (const auto& issue : std::ranges::equal_range(input_warnings, molecule_index, {},
                                                      &calculation::InputIssue::molecule_index)) {
        auto diagnostic = ResultDiagnostic{};
        diagnostic.severity = DiagnosticSeverity::warning;
        diagnostic.code = calculation::to_string(issue.kind);
        diagnostic.message = issue.message;
        diagnostic.molecule_index = issue.molecule_index;
        diagnostic.conformer_index = issue.conformer_index;
        append_unique(diagnostics, seen, std::move(diagnostic));
    }
    if (execution != nullptr) {
        for (const auto& issue :
             std::ranges::equal_range(plan_warnings(*execution), std::optional{molecule_index}, {},
                                      &methods::ExecutionIssue::molecule_index)) {
            append_unique(diagnostics, seen, plan_warning_diagnostic(issue));
        }
    }
    if (const auto diagnostic = status_diagnostic(result); diagnostic.has_value()) {
        append_unique(diagnostics, seen, *diagnostic);
    }
    if (result.status() != ResultStatus::no_executable_plan) {
        return diagnostics;
    }
    auto positions = std::vector<ChargeCalculationResult::IssuePosition>{};
    std::ranges::merge(result.molecule_issues_[molecule_index], result.shared_issues_,
                       std::back_inserter(positions));
    for (const auto& [rejection_index, issue_index] : positions) {
        const auto& rejected = execution->rejections[rejection_index];
        const auto& issue_value = rejected.issues[issue_index];
        auto candidate = "method '" + rejected.method_id + "'";
        if (rejected.parameter_set_id.has_value()) {
            candidate += ", parameter set '" + *rejected.parameter_set_id + "'";
        }
        if (const auto* issue = std::get_if<methods::PrerequisiteIssue>(&issue_value)) {
            auto diagnostic = ResultDiagnostic{};
            diagnostic.severity = DiagnosticSeverity::error;
            diagnostic.code = prerequisite_code(issue->kind);
            diagnostic.message = candidate + ": " + issue->message;
            diagnostic.molecule_index = issue->molecule_index;
            diagnostic.atom_index = issue->atom_index;
            diagnostic.bond_index = issue->bond_index;
            diagnostic.conformer_index = issue->conformer_index;
            append_unique(diagnostics, seen, std::move(diagnostic));
        } else {
            const auto& execution_issue = std::get<methods::ExecutionIssue>(issue_value);
            auto diagnostic = ResultDiagnostic{};
            diagnostic.severity = DiagnosticSeverity::error;
            diagnostic.code = methods::to_string(execution_issue.kind);
            diagnostic.message = candidate + ": " + execution_issue.message;
            diagnostic.molecule_index = execution_issue.molecule_index;
            append_unique(diagnostics, seen, std::move(diagnostic));
        }
    }
    return diagnostics;
}

auto record_assignments(const ChargeCalculationResult& result, const std::size_t record_index)
    -> std::span<const charges::ChargeAssignment> {
    const auto* charge_set = result.charges();
    if (charge_set == nullptr) {
        return {};
    }
    const auto assignments = charge_set->assignments();
    const auto [first, last] =
        std::ranges::equal_range(assignments, record_index, {}, [](const auto& assignment) {
            return assignment.target.molecule_index;
        });
    return {first, last};
}

} // namespace chargefw::adapters
