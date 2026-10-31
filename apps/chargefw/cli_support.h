#pragma once

#include <CLI/CLI.hpp>
#include <chargefw/adapters/gemmi/all.h>
#include <chargefw/adapters/native/all.h>
#include <chargefw/calculation/calculation.h>
#include <chargefw/core/molecule_collection.h>

#include <chrono>
#include <optional>
#include <string>
#include <vector>

namespace chargefw::cli {

struct ImportedCollection {
    core::MoleculeCollection molecules;
    std::vector<adapters::ImportedMoleculeRecord> records;
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
    std::string cutoff_atom_threshold;
    std::string cover_atom_threshold;
    std::size_t max_threads = 0;
    bool permissive_types = false;
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
};

struct CalculationRun {
    adapters::ExecutionMetrics metrics;
    std::chrono::steady_clock::time_point started;
};

[[nodiscard]] auto utc_timestamp() -> std::string;
[[nodiscard]] auto peak_resident_memory_mb() -> double;

void add_input_options(CLI::App& command, InputArguments& arguments);
void add_selection_options(CLI::App& command, SelectionArguments& arguments);
[[nodiscard]] auto import_input(const InputArguments& arguments) -> ImportedCollection;
[[nodiscard]] auto make_request(core::MoleculeCollection molecules,
                                const SelectionArguments& arguments,
                                std::span<const adapters::ImportedMoleculeRecord> records)
    -> calculation::AssessmentRequest;
[[nodiscard]] auto make_requested_provenance(const calculation::AssessmentRequest& request,
                                             std::size_t max_threads)
    -> adapters::RequestedCalculationProvenance;

void print_inspection(const ImportedCollection& imported);
void print_applicability(const calculation::AssessmentResult& assessment, std::size_t max_threads);
void print_methods(const std::string& method_id);
void print_parameter_sets(const std::string& parameter_set_id, const std::string& method_id,
                          bool keys);

[[nodiscard]] auto
write_calculation_outputs(const std::string& output_directory, const std::string& input_path,
                          std::vector<adapters::ImportedMoleculeRecord> records,
                          const adapters::RequestedCalculationProvenance& requested,
                          calculation::ExecutionResult result,
                          const OutputArguments& output_arguments, CalculationRun& run) -> int;

} // namespace chargefw::cli
