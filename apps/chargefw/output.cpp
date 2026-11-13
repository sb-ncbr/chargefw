#include "cli_support.h"
#include "core/diagnostic_description.h"
#include "platform/peak_memory.h"

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <new>
#include <print>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace chargefw::cli {
namespace {

void report_diagnostics(const adapters::ChargeCalculationResult& result,
                        std::set<std::pair<std::string, std::string>>& reported) {
    const auto print = [](const adapters::ResultDiagnostic& diagnostic,
                          const std::string_view message) {
        const auto label = diagnostic.severity == adapters::DiagnosticSeverity::warning ? "Warning"
                           : diagnostic.severity == adapters::DiagnosticSeverity::info  ? "Info"
                                                                                        : "Error";
        std::println(std::cerr, "{}: {}", label, message);
    };
    const auto report = [&](const adapters::ResultDiagnostic& diagnostic) {
        if (reported.emplace(diagnostic.code, diagnostic.message).second) {
            print(diagnostic, diagnostic.message);
        }
    };
    for (const auto& diagnostic : adapters::charge_result_diagnostics(result)) {
        report(diagnostic);
    }

    // Warnings about molecules name the molecule and are grouped by code, so that a file of many
    // records gives one line per kind: the first molecule, or the count with the first example.
    struct WarningGroup {
        adapters::ResultDiagnostic first;
        std::set<std::size_t> molecules;
    };
    auto groups = std::vector<WarningGroup>{};
    auto errors = std::vector<adapters::ResultDiagnostic>{};
    const auto records = result.inputs();
    for (std::size_t index = 0; index < records.size(); ++index) {
        for (auto& diagnostic : adapters::charge_record_diagnostics(result, index)) {
            if (diagnostic.severity != adapters::DiagnosticSeverity::warning ||
                !diagnostic.molecule_index.has_value()) {
                errors.push_back(std::move(diagnostic));
                continue;
            }
            if (reported.contains({diagnostic.code, diagnostic.message})) {
                continue;
            }
            const auto group = std::ranges::find(
                groups, diagnostic.code, [](const auto& value) { return value.first.code; });
            if (group == groups.end()) {
                groups.push_back({.first = diagnostic, .molecules = {index}});
            } else {
                group->molecules.insert(index);
            }
        }
    }
    for (const auto& [first, molecules] : groups) {
        const auto molecule_index = *first.molecule_index;
        const auto example =
            core::detail::molecule_description(records[molecule_index].molecule, molecule_index) +
            ": " + first.message;
        print(first, molecules.size() == 1
                         ? example
                         : std::format("{} of {} molecules: {}, for example {}", molecules.size(),
                                       records.size(), first.code, example));
    }
    for (const auto& diagnostic : errors) {
        report(diagnostic);
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

// Writes beside the destination and replaces it only once writing has succeeded, so a refused or
// failed write leaves an existing file unchanged.
template <typename Writer> void write_output(const std::filesystem::path& path, Writer writer) {
    // A random name, created exclusively: an existing file or symbolic link is never opened.
    auto random = std::random_device{};
    auto temporary = path;
    temporary += std::format(".{:08x}{:08x}.tmp", random(), random());
    auto output = std::ofstream{temporary, std::ios::out | std::ios::noreplace};
    if (!output) {
        throw std::runtime_error{"Unable to open output file: " + path.string()};
    }
    try {
        writer(output);
        finalize_output(output, path);
        std::filesystem::rename(temporary, path);
    } catch (...) {
        // Close first: an open file cannot be removed on Windows.
        output.close();
        auto ignored = std::error_code{};
        std::filesystem::remove(temporary, ignored);
        throw;
    }
}

void write_json(const std::filesystem::path& path, const adapters::ChargeCalculationResult& result,
                const adapters::ExecutionMetrics& metrics) {
    write_output(path, [&result, &metrics](auto& output) {
        adapters::native::json_output::JsonWriter{output}.write(result, metrics);
    });
}

void write_mmcif(const std::filesystem::path& path,
                 const adapters::ChargeCalculationResult& result) {
    write_output(path, [&result](auto& output) {
        adapters::gemmi::mmcif_output::MmcifWriter{output}.write(result);
    });
}

void write_mol2(const std::filesystem::path& path,
                const adapters::ChargeCalculationResult& result) {
    write_output(path, [&result](auto& output) {
        adapters::native::mol2_output::Mol2Writer{output}.write(result);
    });
}

void write_pqr(const std::filesystem::path& path, const adapters::ChargeCalculationResult& result) {
    write_output(path, [&result](auto& output) {
        adapters::native::pqr_output::PqrWriter{output}.write(result);
    });
}

} // namespace

auto utc_timestamp() -> std::string {
    // system_clock is UTC; millisecond precision makes %T print three decimal places.
    return std::format("{:%FT%T}Z", std::chrono::floor<std::chrono::milliseconds>(
                                        std::chrono::system_clock::now()));
}

auto write_calculation_outputs(const std::string& output_directory, const std::string& input_path,
                               const adapters::ChargeCalculationResult& owned_result,
                               const OutputArguments& output_arguments, CalculationRun& run)
    -> int {
    // Every output that cannot be written, the result JSON too, exits with 6. Running out of memory
    // before the result JSON is written exits with 1.
    const auto write_failed = [](const std::exception& error) -> int {
        std::println(std::cerr, "Write error: {}", error.what());
        return 6;
    };
    const auto directory = std::filesystem::path{output_directory};
    const auto prefix =
        directory / (std::filesystem::path{input_path}.stem().string() + ".chargefw");
    run.metrics.peak_resident_memory_bytes = platform::peak_resident_memory_bytes();
    run.metrics.ended_at = utc_timestamp();
    run.metrics.runtime_seconds =
        std::chrono::duration<double>{std::chrono::steady_clock::now() - run.started}.count();
    const auto json_path = std::filesystem::path{prefix.string() + ".json"};
    try {
        std::error_code directory_error;
        std::filesystem::create_directories(directory, directory_error);
        if (directory_error) {
            throw std::runtime_error{"Unable to create output directory: " + directory.string() +
                                     ": " + directory_error.message()};
        }
        if (!std::filesystem::is_directory(directory)) {
            throw std::runtime_error{"Output path is not a directory: " + directory.string()};
        }
        write_json(json_path, owned_result, run.metrics);
    } catch (const std::bad_alloc&) {
        throw;
    } catch (const std::exception& error) {
        return write_failed(error);
    }
    report_diagnostics(owned_result, run.reported);

    if (owned_result.charges() == nullptr) {
        switch (owned_result.status()) {
        case adapters::ResultStatus::invalid_input_or_request:
            return 2;
        case adapters::ResultStatus::no_executable_plan:
            return 3;
        case adapters::ResultStatus::numerical_failure:
            return 4;
        case adapters::ResultStatus::cancelled:
            return 5;
        case adapters::ResultStatus::success:
            throw std::logic_error{"successful result has no charges"};
        }
    }

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
        if (output_arguments.pqr) {
            const auto pqr_path = std::filesystem::path{prefix.string() + ".pqr"};
            write_pqr(pqr_path, owned_result);
            std::println("Wrote {}", pqr_path.string());
        }
    } catch (const std::exception& error) {
        // Also for running out of memory: the result JSON is written, which status 1 rules out.
        return write_failed(error);
    }
    return 0;
}

} // namespace chargefw::cli
