#pragma once

#include <chargefw/adapters/molecule_record.h>
#include <chargefw/calculation/calculation.h>
#include <chargefw/methods/method_options.h>

#include <cstddef>
#include <cstdint>
#include <map>
#include <optional>
#include <span>
#include <string>
#include <vector>

namespace chargefw::adapters {

enum class DiagnosticSeverity : std::uint8_t { info, warning, error };

struct ResultDiagnostic {
    DiagnosticSeverity severity = DiagnosticSeverity::error;
    std::string code;
    std::string message;
    // Structured indices are zero-based. Source line numbers are one-based.
    std::optional<std::size_t> molecule_index;
    std::optional<std::size_t> atom_index;
    std::optional<std::size_t> bond_index;
    std::optional<std::size_t> conformer_index;
    std::optional<std::size_t> line;
};

// Invocation-wide calculation provenance. The JSON writer serializes requested and effective
// policy; molecular output formats do not consume it.
struct RequestedCalculationProvenance {
    std::optional<std::string> method_id;
    std::optional<std::string> parameter_set_id;
    std::string method_scope = "recommended";
    bool permissive_types = false;
    std::optional<std::size_t> cutoff_atom_threshold;
    std::optional<std::size_t> cover_atom_threshold;
    std::size_t max_threads = 0;
    std::string execution_kind = "auto";
    std::optional<double> execution_radius;
    std::map<std::string, methods::MethodOptions> method_options;
};

struct ExecutionMetrics {
    std::string started_at;
    std::string ended_at;
    double runtime_seconds = 0.0;
    double parsing_seconds = 0.0;
    double applicability_seconds = 0.0;
    double computation_seconds = 0.0;
    double peak_resident_memory_mb = 0.0;
};

// Application-facing result boundary. It owns the exact normalized records supplied to the
// calculation together with the requested and effective calculation outcome.
class ChargeCalculationResult {
  public:
    [[nodiscard]] auto inputs() const noexcept -> std::span<const ImportedMoleculeRecord>;
    [[nodiscard]] auto requested() const noexcept -> const RequestedCalculationProvenance&;
    [[nodiscard]] auto execution() const noexcept -> const calculation::ExecutionResult&;

  private:
    // Results come only from make_charge_calculation_result(), which validates them.
    ChargeCalculationResult() = default;

    friend auto make_charge_calculation_result(std::vector<ImportedMoleculeRecord> inputs,
                                               RequestedCalculationProvenance requested,
                                               calculation::ExecutionResult execution)
        -> ChargeCalculationResult;

    std::vector<ImportedMoleculeRecord> inputs_;
    RequestedCalculationProvenance requested_;
    calculation::ExecutionResult execution_;
};

[[nodiscard]] auto make_charge_calculation_result(std::vector<ImportedMoleculeRecord> inputs,
                                                  RequestedCalculationProvenance requested,
                                                  calculation::ExecutionResult execution)
    -> ChargeCalculationResult;

[[nodiscard]] auto charge_result_diagnostics(const ChargeCalculationResult& result)
    -> std::vector<ResultDiagnostic>;

[[nodiscard]] auto charge_record_diagnostics(const ChargeCalculationResult& result,
                                             std::size_t molecule_index)
    -> std::vector<ResultDiagnostic>;

} // namespace chargefw::adapters
