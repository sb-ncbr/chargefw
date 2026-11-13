#pragma once

#include <chargefw/calculation/calculation.h>
#include <chargefw/charges/charge_set.h>
#include <chargefw/core/issue_location.h>
#include <chargefw/io/molecule_record.h>
#include <chargefw/methods/method_options.h>

#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

namespace chargefw::io {

enum class DiagnosticSeverity : std::uint8_t {
    info,
    warning,
    error,
};

// The outcome an application reports for its request. Unlike calculation::ExecutionStatus, it
// includes a request that failed with an exception, so that no execution result exists.
enum class ResultStatus : std::uint8_t {
    success,
    invalid_input_or_request,
    no_executable_plan,
    numerical_failure,
    cancelled,
};

// The value result JSON reports, such as "invalid_input_or_request".
[[nodiscard]] auto to_string(ResultStatus value) -> std::string_view;

struct ResultDiagnostic {
    DiagnosticSeverity severity = DiagnosticSeverity::error;
    std::string code;
    std::string message;
    // Structured indices are zero-based. Source line numbers are one-based.
    std::optional<std::size_t> molecule_index;
    core::IssueLocation location;
    std::optional<std::size_t> line;
};

// Invocation-wide calculation provenance. The JSON writer serializes requested and effective
// policy; molecular output formats do not consume it.
struct RequestedCalculationProvenance {
    std::optional<std::string> method_id;
    std::optional<std::string> parameter_set_id;
    calculation::MethodScope method_scope = calculation::MethodScope::recommended;
    parameters::ParameterMatching parameter_matching = parameters::ParameterMatching::strict;
    std::optional<std::size_t> cutoff_threshold;
    std::optional<std::size_t> cover_threshold;
    std::size_t threads = 0;
    calculation::ExecutionSelectionMode execution = calculation::ExecutionSelectionMode::automatic;
    std::optional<double> radius;
    // The options given for the requested method.
    methods::MethodOptions method_options;
    // The residue names requested as fixed ions, also those not found in the molecules.
    std::vector<std::string> fixed_ions;
};

struct ExecutionMetrics {
    std::string started_at;
    std::string ended_at;
    double runtime_seconds = 0.0;
    double parsing_seconds = 0.0;
    double applicability_seconds = 0.0;
    double computation_seconds = 0.0;
    std::uint64_t peak_resident_memory_bytes = 0;
};

// Application-facing result boundary. It owns the exact normalized records supplied to the
// calculation together with the requested and effective calculation outcome.
class ChargeCalculationResult {
  public:
    [[nodiscard]] auto inputs() const noexcept -> std::span<const ImportedMoleculeRecord>;
    [[nodiscard]] auto requested() const noexcept -> const RequestedCalculationProvenance&;
    // The calculation's status, or invalid_input_or_request when the request failed.
    [[nodiscard]] auto status() const noexcept -> ResultStatus;
    // Null when the request failed with an exception, so that no execution result exists.
    [[nodiscard]] auto execution() const noexcept -> const calculation::ExecutionResult*;
    // The message of the exception with which the request failed.
    [[nodiscard]] auto request_failure() const noexcept -> std::optional<std::string_view>;
    // The calculated charges; null unless the calculation succeeded.
    [[nodiscard]] auto charges() const noexcept -> const charges::ChargeSet*;

  private:
    struct RequestFailure {
        std::string message;
    };

    // Results come only from the two factories below, each of which creates one alternative.
    ChargeCalculationResult() = default;

    friend auto make_charge_calculation_result(std::vector<ImportedMoleculeRecord> inputs,
                                               RequestedCalculationProvenance requested,
                                               calculation::ExecutionResult execution)
        -> ChargeCalculationResult;
    friend auto make_failed_request_result(std::vector<ImportedMoleculeRecord> inputs,
                                           RequestedCalculationProvenance requested,
                                           std::string message) -> ChargeCalculationResult;
    friend auto charge_record_diagnostics(const ChargeCalculationResult& result,
                                          std::size_t molecule_index)
        -> std::vector<ResultDiagnostic>;

    // Positions (rejection index, issue index) of rejection issues in execution order, so that
    // record diagnostics need not scan the issues of every other molecule. Built only when no
    // plan was executable, the only case in which records report rejections.
    using IssuePosition = std::pair<std::size_t, std::size_t>;

    std::vector<ImportedMoleculeRecord> inputs_;
    RequestedCalculationProvenance requested_;
    std::variant<calculation::ExecutionResult, RequestFailure> outcome_;
    std::vector<std::vector<IssuePosition>> molecule_issues_;
    std::vector<IssuePosition> shared_issues_;
};

[[nodiscard]] auto make_charge_calculation_result(std::vector<ImportedMoleculeRecord> inputs,
                                                  RequestedCalculationProvenance requested,
                                                  calculation::ExecutionResult execution)
    -> ChargeCalculationResult;

// What was requested, for result JSON: the choices of `request`, the thread limit given to
// calculate(), and the residue names the fixed ions were selected by.
[[nodiscard]] auto make_requested_provenance(const calculation::AssessmentRequest& request,
                                             std::size_t threads,
                                             std::vector<std::string> fixed_ions = {})
    -> RequestedCalculationProvenance;

// For a request that failed with an exception after its records and options were available.
[[nodiscard]] auto make_failed_request_result(std::vector<ImportedMoleculeRecord> inputs,
                                              RequestedCalculationProvenance requested,
                                              std::string message) -> ChargeCalculationResult;

[[nodiscard]] auto charge_result_diagnostics(const ChargeCalculationResult& result)
    -> std::vector<ResultDiagnostic>;

[[nodiscard]] auto charge_record_diagnostics(const ChargeCalculationResult& result,
                                             std::size_t molecule_index)
    -> std::vector<ResultDiagnostic>;

} // namespace chargefw::io
