#include "cli_support.h"

#include "adapters/fixed_ions.h"

#include <chargefw/methods/method_registry.h>
#include <chargefw/parameters/io/parameter_set_io.h>

#include <algorithm>
#include <cctype>
#include <charconv>
#include <filesystem>
#include <fstream>
#include <optional>
#include <ranges>
#include <stdexcept>
#include <utility>
#include <variant>

namespace chargefw::cli {
namespace {

template <typename Reader>
auto read_collection(Reader& reader, const std::string& input_path) -> ImportedCollection {
    std::vector<core::Molecule> molecules;
    std::vector<adapters::ImportedMoleculeRecord> records;

    while (auto record = reader.next()) {
        auto imported = std::move(*record);
        molecules.push_back(imported.molecule);
        records.push_back(std::move(imported));
    }

    if (molecules.empty()) {
        throw std::runtime_error{"No molecules found in input file: " + input_path};
    }

    return ImportedCollection{.molecules =
                                  core::MoleculeCollection{std::move(molecules), input_path},
                              .records = std::move(records)};
}

[[nodiscard]] auto parse_atom_threshold(const std::string& value, const std::string_view name)
    -> std::optional<std::size_t> {
    if (value == "unlimited") {
        return std::nullopt;
    }

    std::size_t threshold = 0;
    const auto [end, error] = std::from_chars(value.data(), value.data() + value.size(), threshold);
    if (error != std::errc{} || end != value.data() + value.size()) {
        throw std::runtime_error{std::string{name} +
                                 " must be a non-negative integer or 'unlimited'"};
    }
    return threshold;
}

[[nodiscard]] auto parse_method_option(const std::string& text)
    -> std::pair<std::string, std::pair<std::string, methods::MethodOptionValue>> {
    const auto equals = text.find('=');
    const auto dot = text.find('.');
    if (dot == std::string::npos || equals == std::string::npos || dot == 0 || dot + 1 >= equals ||
        equals + 1 >= text.size()) {
        throw std::runtime_error{"Method option must have the form METHOD.OPTION=VALUE: " + text};
    }

    const auto method_id = text.substr(0, dot);
    const auto option_id = text.substr(dot + 1, equals - dot - 1);
    const auto value = text.substr(equals + 1);
    const auto* method = methods::method_registry().find(method_id);
    if (method == nullptr) {
        throw std::runtime_error{"Unknown method in method option: " + method_id};
    }
    const auto schema = method->option_schema();
    const auto spec =
        std::ranges::find_if(schema, [&option_id](const methods::MethodOptionSpec& item) {
            return item.id == option_id;
        });
    if (spec == schema.end()) {
        throw std::runtime_error{"Unknown option '" + option_id + "' for method '" + method_id +
                                 "'"};
    }

    switch (spec->type) {
    case methods::MethodOptionType::boolean:
        if (value == "true" || value == "1") {
            return {method_id, {option_id, true}};
        }
        if (value == "false" || value == "0") {
            return {method_id, {option_id, false}};
        }
        break;
    case methods::MethodOptionType::integer: {
        int parsed = 0;
        const auto [end, error] =
            std::from_chars(value.data(), value.data() + value.size(), parsed);
        if (error == std::errc{} && end == value.data() + value.size()) {
            return {method_id, {option_id, parsed}};
        }
        break;
    }
    case methods::MethodOptionType::floating_point: {
        try {
            std::size_t end = 0;
            const auto parsed = std::stod(value, &end);
            if (end == value.size()) {
                return {method_id, {option_id, parsed}};
            }
        } catch (const std::exception&) {
            break;
        }
        break;
    }
    case methods::MethodOptionType::string:
        return {method_id, {option_id, value}};
    }
    throw std::runtime_error{"Invalid value for method option '" + method_id + "." + option_id +
                             "': " + value};
}

auto read_collection(const std::string& input_path,
                     const adapters::gemmi::InputOptions structural_options,
                     const bool structural_options_requested) -> ImportedCollection {
    std::ifstream input{input_path};
    if (!input) {
        throw std::runtime_error{"Unable to open input file: " + input_path};
    }

    auto extension = std::filesystem::path{input_path}.extension().string();
    std::ranges::transform(extension, extension.begin(), [](const unsigned char character) -> char {
        return static_cast<char>(std::tolower(character));
    });
    if (extension == ".sdf") {
        if (structural_options_requested) {
            throw std::runtime_error{
                "Structural input options are only supported for PDB and mmCIF input"};
        }
        adapters::native::sdf_input::SdfReader reader{input, input_path};
        return read_collection(reader, input_path);
    }
    if (extension == ".mol") {
        if (structural_options_requested) {
            throw std::runtime_error{
                "Structural input options are only supported for PDB and mmCIF input"};
        }
        adapters::native::mol_input::MolReader reader{input, input_path};
        return read_collection(reader, input_path);
    }
    if (extension == ".mol2") {
        if (structural_options_requested) {
            throw std::runtime_error{
                "Structural input options are only supported for PDB and mmCIF input"};
        }
        adapters::native::mol2_input::Mol2Reader reader{input, input_path};
        return read_collection(reader, input_path);
    }
    if (extension == ".json") {
        if (structural_options_requested) {
            throw std::runtime_error{
                "Structural input options are only supported for PDB and mmCIF input"};
        }
        adapters::native::json_input::JsonReader reader{input, input_path,
                                                        structural_options.conformers};
        return read_collection(reader, input_path);
    }
    if (extension == ".pdb") {
        adapters::gemmi::pdb_input::PdbReader reader{input, input_path, structural_options};
        return read_collection(reader, input_path);
    }
    if (extension == ".cif" || extension == ".mmcif") {
        adapters::gemmi::mmcif_input::MmcifReader reader{input, input_path, structural_options};
        return read_collection(reader, input_path);
    }

    throw std::runtime_error{"Unsupported input file type: " + extension +
                             ". Supported types: .sdf, .mol, .mol2, .json, .pdb, .cif, .mmcif"};
}

} // namespace

void add_input_options(CLI::App& command, InputArguments& arguments) {
    command.add_option("input", arguments.path, "Input molecular file")->required();
    arguments.structural_selection_option =
        command.add_option("--structural-selection", arguments.structural_selection,
                           "PDB/mmCIF record selection: all, polymers-and-ligands, or polymers");
    arguments.structural_bonds_option =
        command.add_option("--structural-bonds", arguments.structural_bonds,
                           "PDB/mmCIF connectivity: none, explicit, templates, or hybrid");
    command.add_option("--conformers", arguments.conformer_selection,
                       "Conformers to read: first or all");
}

void add_selection_options(CLI::App& command, SelectionArguments& arguments) {
    arguments.method_option = command.add_option("--method", arguments.method_id, "Method ID");
    arguments.parameter_set_option =
        command.add_option("--parameter-set", arguments.parameter_set_id, "Parameter-set ID");
    arguments.parameter_set_option->needs(arguments.method_option);
    command
        .add_option("--method-scope", arguments.method_scope,
                    "Methods automatic selection may use: recommended or all")
        ->excludes(arguments.method_option);
    command.add_option("--fixed-ions", arguments.fixed_ion_ids,
                       "Fixed-ion component ID (repeatable)");
    command.add_flag("--permissive-types", arguments.permissive_types,
                     "Allow permissive parameter type classification");
    command
        .add_option("--method-option", arguments.method_options,
                    "Method option METHOD.OPTION=VALUE (repeatable)")
        ->allow_extra_args(false);
    command.add_option("--execution", arguments.execution,
                       "Execution: auto, full, cutoff, or cover");
    command.add_option("--radius", arguments.radius, "Cutoff or cover radius in angstrom");
    arguments.cutoff_atom_threshold_option =
        command.add_option("--cutoff-atom-threshold", arguments.cutoff_atom_threshold,
                           "Automatic full-to-cutoff atom threshold, or unlimited");
    arguments.cover_atom_threshold_option =
        command.add_option("--cover-atom-threshold", arguments.cover_atom_threshold,
                           "Automatic cutoff-to-cover atom threshold, or unlimited");
}

auto import_input(const InputArguments& arguments) -> ImportedCollection {
    const auto options = adapters::gemmi::InputOptions{
        .selection = adapters::gemmi::record_selection_from_string(arguments.structural_selection),
        .bond_strategy = adapters::gemmi::bond_strategy_from_string(arguments.structural_bonds),
        .conformers = adapters::conformer_selection_from_string(arguments.conformer_selection)};
    auto imported = read_collection(arguments.path, options,
                                    arguments.structural_selection_option->count() > 0 ||
                                        arguments.structural_bonds_option->count() > 0);
    return imported;
}

auto make_request(core::MoleculeCollection molecules, const SelectionArguments& arguments,
                  const std::span<const adapters::ImportedMoleculeRecord> records)
    -> calculation::AssessmentRequest {
    auto method_options = std::unordered_map<std::string, methods::MethodOptions>{};
    for (const auto& text : arguments.method_options) {
        const auto [method_id, option] = parse_method_option(text);
        method_options[method_id].set(option.first, option.second);
    }
    auto request = calculation::AssessmentRequest{
        .molecules = std::move(molecules),
        .parameter_sets = parameters::load_default_parameter_sets(),
        .method_id = arguments.method_option->count() == 0 ? std::nullopt
                                                           : std::optional{arguments.method_id},
        .parameter_set_id = arguments.parameter_set_option->count() == 0
                                ? std::nullopt
                                : std::optional{arguments.parameter_set_id},
        .method_scope = calculation::method_scope_from_string(arguments.method_scope),
        .method_options = std::move(method_options),
        .classification_options = {.permissive_types = arguments.permissive_types},
        .execution_selection =
            calculation::ExecutionSelection{
                calculation::execution_selection_kind_from_string(arguments.execution),
                arguments.radius},
        .resource_policy = {
            .cutoff_atom_threshold =
                arguments.cutoff_atom_threshold_option->count() == 0
                    ? std::optional<std::size_t>{calculation::default_cutoff_atom_threshold}
                    : parse_atom_threshold(arguments.cutoff_atom_threshold,
                                           "Cutoff atom threshold"),
            .cover_atom_threshold =
                arguments.cover_atom_threshold_option->count() == 0
                    ? std::optional<std::size_t>{calculation::default_cover_atom_threshold}
                    : parse_atom_threshold(arguments.cover_atom_threshold,
                                           "Cover atom threshold")}};
    if (!arguments.fixed_ion_ids.empty()) {
        auto resolved = adapters::detail::resolve_fixed_ions(records, arguments.fixed_ion_ids);
        if (!resolved.sources.empty()) {
            request.fixed_ions = std::move(resolved);
        }
    }
    return request;
}

} // namespace chargefw::cli
