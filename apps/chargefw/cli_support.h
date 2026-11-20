#pragma once

#include <CLI/CLI.hpp>
#include <chargefw/calculation/calculation.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/io/all.h>

#include <chrono>
#include <cstddef>
#include <filesystem>
#include <optional>
#include <set>
#include <span>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace chargefw::cli {

struct ImportedCollection {
    core::MoleculeCollection molecules;
    std::vector<io::ImportedMoleculeRecord> records;
};

struct InputArguments {
    std::string path;
    std::string structural_selection = "all";
    std::string structural_bonds = "hybrid";
    std::string conformer_selection = "all";
    CLI::Option* structural_selection_option = nullptr;
    CLI::Option* structural_bonds_option = nullptr;
};

struct SelectionArguments {
    std::string method_id;
    std::string parameter_set_id;
    std::string method_scope = "recommended";
    std::string execution = "auto";
    std::optional<double> radius;
    std::string cutoff_threshold;
    std::string cover_threshold;
    std::size_t max_threads = 0;
    std::string parameter_matching = "strict";
    std::vector<std::string> method_options;
    std::vector<std::string> fixed_ion_ids;
    CLI::Option* method_option = nullptr;
    CLI::Option* parameter_set_option = nullptr;
    CLI::Option* cutoff_atom_threshold_option = nullptr;
    CLI::Option* cover_atom_threshold_option = nullptr;
};

struct OutputArguments {
    bool mol2 = false;
    bool mmcif = false;
    bool pqr = false;
    bool annotate_mmcif = false;
    std::string if_exists = "error";
};

struct CalculationRun {
    io::ExecutionMetrics metrics;
    std::chrono::steady_clock::time_point started;
    // Diagnostics already printed, by code and message, so that they are not printed again.
    std::set<std::pair<std::string, std::string>> reported;
};

[[nodiscard]] auto utc_timestamp() -> std::string;

// One line for the warnings of one code about molecules: the molecule when there is one, or the
// count with the first one as an example. `molecules` is the zero-based indices, in input order;
// `conformer_index` is that of the first one's warning.
[[nodiscard]] auto molecule_warning_text(std::span<const io::ImportedMoleculeRecord> records,
                                         std::string_view code, std::string_view message,
                                         const std::set<std::size_t>& molecules,
                                         std::optional<std::size_t> conformer_index = std::nullopt)
    -> std::string;

// The input path without a final .gz extension, in any case, which names a gzip-compressed file.
[[nodiscard]] auto uncompressed_path(const std::string& path) -> std::filesystem::path;

void add_input_options(CLI::App& command, InputArguments& arguments);
void add_selection_options(CLI::App& command, SelectionArguments& arguments);
[[nodiscard]] auto import_input(const InputArguments& arguments) -> ImportedCollection;
[[nodiscard]] auto make_request(core::MoleculeCollection molecules,
                                const SelectionArguments& arguments,
                                std::span<const io::ImportedMoleculeRecord> records)
    -> calculation::AssessmentRequest;

void print_inspection(const ImportedCollection& imported);
void print_applicability(const calculation::AssessmentResult& assessment, std::size_t max_threads);
void print_methods(const std::string& method_id);
void print_parameter_sets(const std::string& parameter_set_id, const std::string& method_id,
                          bool keys);

// Refuses, before the calculation, records that --annotate-mmcif cannot write.
void validate_annotation(std::span<const io::ImportedMoleculeRecord> records,
                         const OutputArguments& arguments);
[[nodiscard]] auto write_calculation_outputs(const std::string& output_directory,
                                             const std::string& input_path,
                                             const io::ChargeCalculationResult& owned_result,
                                             const OutputArguments& output_arguments,
                                             CalculationRun& run) -> int;

} // namespace chargefw::cli
