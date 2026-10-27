#include "cli_support.h"
#include "platform/peak_memory.h"

#include <chargefw/config.h>

#include <chrono>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <print>
#include <set>
#include <stdexcept>

namespace chargefw::cli {
namespace {

void report_diagnostics(const adapters::ChargeCalculationResult& result) {
    auto reported = std::set<std::pair<std::string, std::string>>{};
    const auto report = [&reported](const adapters::ResultDiagnostic& diagnostic) {
        if (!reported.emplace(diagnostic.code, diagnostic.message).second) {
            return;
        }
        const auto label = diagnostic.severity == adapters::DiagnosticSeverity::warning ? "Warning"
                           : diagnostic.severity == adapters::DiagnosticSeverity::info  ? "Info"
                                                                                        : "Error";
        std::println(std::cerr, "{}: {}", label, diagnostic.message);
    };
    for (const auto& diagnostic : adapters::charge_result_diagnostics(result)) {
        report(diagnostic);
    }
    for (std::size_t index = 0; index < result.inputs().size(); ++index) {
        for (const auto& diagnostic : adapters::charge_record_diagnostics(result, index)) {
            report(diagnostic);
        }
    }
}

void finalize_output(std::ofstream& output, const std::filesystem::path& path) {
    if (!output) {
        throw std::runtime_error{"Unable to write output file: " + path.string()};
    }
    output.flush();
    if (!output) {
        throw std::runtime_error{"Unable to write output file: " + path.string()};
    }
    output.close();
    if (!output) {
        throw std::runtime_error{"Unable to write output file: " + path.string()};
    }
}

template <typename Writer> void write_output(const std::filesystem::path& path, Writer writer) {
    auto output = std::ofstream{path};
    if (!output) {
        throw std::runtime_error{"Unable to open output file: " + path.string()};
    }
    writer(output);
    finalize_output(output, path);
}

void write_json(const std::filesystem::path& path, const adapters::ChargeCalculationResult& result,
                const adapters::ExecutionMetrics& metrics) {
    write_output(path, [&result, &metrics](auto& output) {
        adapters::native::json_output::JsonWriter{output}.write(result, CHARGEFW_VERSION_STRING,
                                                                metrics);
    });
}

void write_mmcif(const std::filesystem::path& path,
                 const adapters::ChargeCalculationResult& result) {
    write_output(path, [&result](auto& output) {
        adapters::gemmi::mmcif_output::MmcifWriter{output}.write(result, CHARGEFW_VERSION_STRING);
    });
}

void write_mol2(const std::filesystem::path& path,
                const adapters::ChargeCalculationResult& result) {
    write_output(path, [&result](auto& output) {
        adapters::native::mol2_output::Mol2Writer{output}.write(result, CHARGEFW_VERSION_STRING);
    });
}

} // namespace

auto utc_timestamp() -> std::string {
    // system_clock is UTC; millisecond precision makes %T print three decimal places.
    return std::format("{:%FT%T}Z", std::chrono::floor<std::chrono::milliseconds>(
                                        std::chrono::system_clock::now()));
}

auto peak_resident_memory_mb() -> double {
    return static_cast<double>(platform::peak_resident_memory_bytes()) / (1024.0 * 1024.0);
}

auto make_requested_provenance(const calculation::AssessmentRequest& request,
                               const std::size_t max_threads)
    -> adapters::RequestedCalculationProvenance {
    auto requested = adapters::RequestedCalculationProvenance{
        .method_id = request.method_id,
        .parameter_set_id = request.parameter_set_id,
        .permissive_types = request.classification_options.permissive_types,
        .cutoff_atom_threshold = request.resource_policy.cutoff_atom_threshold,
        .cover_atom_threshold = request.resource_policy.cover_atom_threshold,
        .max_threads = max_threads,
        .execution_kind = std::string{calculation::to_string(request.execution_selection.kind())},
        .execution_radius = request.execution_selection.radius()};
    for (const auto& [method_id, options] : request.method_options) {
        requested.method_options.emplace(method_id, options);
    }
    return requested;
}

auto write_calculation_outputs(const std::string& output_directory, const std::string& input_path,
                               const std::vector<adapters::ImportedMoleculeRecord>& records,
                               const adapters::RequestedCalculationProvenance& requested,
                               const calculation::ExecutionResult& result,
                               const OutputArguments& output_arguments, CalculationRun& run)
    -> int {
    const auto owned_result = adapters::make_charge_calculation_result(records, requested, result);
    run.metrics.peak_resident_memory_mb = peak_resident_memory_mb();
    const auto directory = std::filesystem::path{output_directory};
    std::error_code directory_error;
    std::filesystem::create_directories(directory, directory_error);
    if (directory_error) {
        throw std::runtime_error{"Unable to create output directory: " + directory.string() + ": " +
                                 directory_error.message()};
    }
    if (!std::filesystem::is_directory(directory)) {
        throw std::runtime_error{"Output path is not a directory: " + directory.string()};
    }
    const auto prefix =
        directory / (std::filesystem::path{input_path}.stem().string() + ".chargefw");
    if (!result.calculated()) {
        run.metrics.ended_at = utc_timestamp();
        run.metrics.runtime_seconds =
            std::chrono::duration<double>{std::chrono::steady_clock::now() - run.started}.count();
        write_json(prefix.string() + ".json", owned_result, run.metrics);
        report_diagnostics(owned_result);
        switch (result.status) {
        case calculation::ExecutionStatus::invalid_input_or_request:
            return 2;
        case calculation::ExecutionStatus::no_executable_plan:
            return 3;
        case calculation::ExecutionStatus::numerical_failure:
            return 4;
        case calculation::ExecutionStatus::cancelled:
            return 5;
        case calculation::ExecutionStatus::success:
            throw std::logic_error{"successful result has no charges"};
        }
    }

    if (!result.charges.has_value()) {
        throw std::runtime_error{"calculation result is missing charges"};
    }
    run.metrics.peak_resident_memory_mb = peak_resident_memory_mb();
    run.metrics.ended_at = utc_timestamp();
    run.metrics.runtime_seconds =
        std::chrono::duration<double>{std::chrono::steady_clock::now() - run.started}.count();
    const auto json_path = std::filesystem::path{prefix.string() + ".json"};
    write_json(json_path, owned_result, run.metrics);
    report_diagnostics(owned_result);
    std::println("Wrote {}", json_path.string());

    try {
        if (output_arguments.mol2) {
            const auto mol2_path = std::filesystem::path{prefix.string() + ".mol2"};
            write_mol2(mol2_path, owned_result);
            std::println("Wrote {}", mol2_path.string());
        }
        if (output_arguments.mmcif) {
            const auto mmcif_path = std::filesystem::path{prefix.string() + ".cif"};
            write_mmcif(mmcif_path, owned_result);
            std::println("Wrote {}", mmcif_path.string());
        }
    } catch (const std::exception& error) {
        std::println(std::cerr, "Export error: {}", error.what());
        return 6;
    }
    return 0;
}

} // namespace chargefw::cli
