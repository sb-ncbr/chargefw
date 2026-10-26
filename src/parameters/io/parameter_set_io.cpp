#include <chargefw/core/periodic_table.h>
#include <chargefw/parameters/io/parameter_data_paths.h>
#include <chargefw/parameters/io/parameter_set_io.h>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace chargefw::parameters {
namespace {

using Json = nlohmann::json;

[[nodiscard]] auto json_type_name(const Json& value) -> std::string {
    return std::string{value.type_name()};
}

[[noreturn]] auto throw_error(const std::string& context, const std::string& message) -> void {
    if (context.empty()) {
        throw std::invalid_argument{message};
    }

    throw std::invalid_argument{context + ": " + message};
}

[[nodiscard]] auto child_context(const std::string& context, const std::string_view child)
    -> std::string {
    if (context.empty()) {
        return std::string{child};
    }

    return context + "." + std::string{child};
}

[[nodiscard]] auto array_context(const std::string& context, const std::size_t index)
    -> std::string {
    return context + "[" + std::to_string(index) + "]";
}

auto ensure_object(const Json& value, const std::string& context) -> void {
    if (!value.is_object()) {
        throw_error(context, "expected object, got " + json_type_name(value));
    }
}

auto ensure_array(const Json& value, const std::string& context) -> void {
    if (!value.is_array()) {
        throw_error(context, "expected array, got " + json_type_name(value));
    }
}

[[nodiscard]] auto required_member(const Json& object, const std::string_view name,
                                   const std::string& context) -> const Json& {
    ensure_object(object, context);

    const auto found = object.find(std::string{name});

    if (found == object.end()) {
        throw_error(context, "missing required member '" + std::string{name} + "'");
    }

    return *found;
}

[[nodiscard]] auto optional_member(const Json& object, const std::string_view name) -> const Json* {
    if (!object.is_object()) {
        return nullptr;
    }

    const auto found = object.find(std::string{name});

    if (found == object.end()) {
        return nullptr;
    }

    return std::addressof(*found);
}

[[nodiscard]] auto require_string(const Json& value, const std::string& context) -> std::string {
    if (!value.is_string()) {
        throw_error(context, "expected string, got " + json_type_name(value));
    }

    return value.get<std::string>();
}

[[nodiscard]] auto require_double(const Json& value, const std::string& context) -> double {
    if (!value.is_number()) {
        throw_error(context, "expected number, got " + json_type_name(value));
    }

    const auto number = value.get<double>();

    if (!std::isfinite(number)) {
        throw_error(context, "expected finite number");
    }

    return number;
}

[[nodiscard]] auto require_uint16(const Json& value, const std::string& context) -> std::uint16_t {
    if (!value.is_number_integer()) {
        throw_error(context, "expected integer, got " + json_type_name(value));
    }

    if (value.is_number_unsigned()) {
        const auto number = value.get<std::uint64_t>();
        if (number > std::numeric_limits<std::uint16_t>::max()) {
            throw_error(context, "expected integer in range [0, 65535]");
        }
        return static_cast<std::uint16_t>(number);
    }

    const auto number = value.get<std::int64_t>();
    if (!std::in_range<std::uint16_t>(number)) {
        throw_error(context, "expected integer in range [0, 65535]");
    }

    return static_cast<std::uint16_t>(number);
}

[[nodiscard]] auto optional_string_member(const Json& object, const std::string_view name,
                                          const std::string& context) -> std::string {
    const auto* value = optional_member(object, name);

    if (value == nullptr) {
        return {};
    }

    return require_string(*value, child_context(context, name));
}

[[nodiscard]] auto parse_metadata(const Json& root, const std::string& context)
    -> ParameterSetMetadata {
    const auto metadata_context = child_context(context, "metadata");
    const auto& metadata_json = required_member(root, "metadata", context);

    ensure_object(metadata_json, metadata_context);

    auto metadata = ParameterSetMetadata{};

    metadata.id = require_string(required_member(metadata_json, "id", metadata_context),
                                 child_context(metadata_context, "id"));

    if (metadata.id.empty()) {
        throw_error(child_context(metadata_context, "id"), "must not be empty");
    }

    metadata.method_id = require_string(required_member(metadata_json, "method", metadata_context),
                                        child_context(metadata_context, "method"));

    if (metadata.method_id.empty()) {
        throw_error(child_context(metadata_context, "method"), "must not be empty");
    }

    metadata.name = optional_string_member(metadata_json, "name", metadata_context);
    metadata.publication = optional_string_member(metadata_json, "publication", metadata_context);
    metadata.notes = optional_string_member(metadata_json, "notes", metadata_context);

    if (const auto* priority = optional_member(metadata_json, "priority")) {
        metadata.priority = require_uint16(*priority, child_context(metadata_context, "priority"));
    }

    return metadata;
}

[[nodiscard]] auto parse_names(const Json& object, const std::string& context)
    -> std::vector<std::string> {
    const auto names_context = child_context(context, "names");
    const auto& names_json = required_member(object, "names", context);

    ensure_array(names_json, names_context);

    std::vector<std::string> names;
    names.reserve(names_json.size());

    for (std::size_t index = 0; index < names_json.size(); ++index) {
        auto name = require_string(names_json[index], array_context(names_context, index));
        if (std::ranges::find(names, name) != names.end()) {
            throw_error(array_context(names_context, index), "duplicate name '" + name + "'");
        }
        names.push_back(std::move(name));
    }

    return names;
}

[[nodiscard]] auto make_named_parameters(const std::vector<std::string>& names,
                                         const Json& values_json, const std::string& context)
    -> std::vector<NamedParameter> {
    ensure_array(values_json, context);

    if (values_json.size() != names.size()) {
        throw_error(context, "value count " + std::to_string(values_json.size()) +
                                 " does not match name count " + std::to_string(names.size()));
    }

    std::vector<NamedParameter> parameters;
    parameters.reserve(names.size());

    for (std::size_t index = 0; index < names.size(); ++index) {
        parameters.push_back(NamedParameter{
            .name = names[index],
            .value = require_double(values_json[index], array_context(context, index))});
    }

    return parameters;
}

[[nodiscard]] auto parse_common_parameters(const Json& root, const std::string& context)
    -> CommonParameters {
    const auto* common_json = optional_member(root, "common");

    if (common_json == nullptr) {
        return CommonParameters{};
    }

    const auto common_context = child_context(context, "common");
    ensure_object(*common_json, common_context);

    const auto names = parse_names(*common_json, common_context);

    const auto values_context = child_context(common_context, "values");
    const auto& values_json = required_member(*common_json, "values", common_context);

    return CommonParameters{make_named_parameters(names, values_json, values_context)};
}

[[nodiscard]] auto require_named_string(const Json& object, const std::string_view name,
                                        const std::string& context) -> std::string {
    const auto field_context = child_context(context, name);
    return require_string(required_member(object, name, context), field_context);
}

[[nodiscard]] auto atomic_number_from_symbol(const std::string& symbol, const std::string& context)
    -> int {
    if (symbol == "*") {
        return 0;
    }

    try {
        return core::periodic_table().element(symbol).atomic_number;
    } catch (const std::out_of_range&) {
        throw_error(context, "unknown element symbol '" + symbol + "'");
    }
}

// Converts a classifier name, reporting an unknown name at its JSON location.
template <typename Convert>
[[nodiscard]] auto classifier_kind(const Convert& convert, const std::string& classifier,
                                   const std::string& context) {
    try {
        return convert(classifier);
    } catch (const std::invalid_argument& error) {
        throw_error(context, error.what());
    }
}

// A plain classifier always derives type "*", so any other type could never match.
auto ensure_plain_type(const bool plain, const std::string& type, const std::string& context)
    -> void {
    if (plain && type != "*") {
        throw_error(child_context(context, "type"), "plain classifier requires type '*'");
    }
}

[[nodiscard]] auto parse_atom_key(const Json& key_json, const std::string& context)
    -> AtomParameterKey {
    ensure_object(key_json, context);

    const auto element = require_named_string(key_json, "element", context);
    const auto classifier = require_named_string(key_json, "classifier", context);
    const auto type = require_named_string(key_json, "type", context);

    const auto classification = classifier_kind(atom_classification_kind_from_string, classifier,
                                                child_context(context, "classifier"));
    ensure_plain_type(classification == AtomParameterClassificationKind::PLAIN, type, context);

    return AtomParameterKey{
        .atomic_number = atomic_number_from_symbol(element, child_context(context, "element")),
        .classification = classification,
        .type = type};
}

[[nodiscard]] auto parse_bond_type_key(const Json& key_json, const std::string& context)
    -> BondTypeKey {
    ensure_object(key_json, context);

    const auto classifier = require_named_string(key_json, "classifier", context);
    const auto type = require_named_string(key_json, "type", context);

    const auto classification = classifier_kind(bond_classification_kind_from_string, classifier,
                                                child_context(context, "classifier"));
    ensure_plain_type(classification == BondParameterClassificationKind::PLAIN, type, context);

    return BondTypeKey{.classification = classification, .type = type};
}

[[nodiscard]] auto parse_bond_key(const Json& key_json, const std::string& context)
    -> BondParameterKey {
    const auto atoms_context = child_context(context, "atoms");
    const auto& atoms_json = required_member(key_json, "atoms", context);
    ensure_array(atoms_json, atoms_context);

    if (atoms_json.size() != 2) {
        throw_error(atoms_context, "must have exactly 2 atom objects");
    }

    const auto bond_context = child_context(context, "bond");
    return BondParameterKey{
        .first_atom = parse_atom_key(atoms_json[0], array_context(atoms_context, 0)),
        .second_atom = parse_atom_key(atoms_json[1], array_context(atoms_context, 1)),
        .bond = parse_bond_type_key(required_member(key_json, "bond", context), bond_context)};
}

[[nodiscard]] auto parse_atom_parameters(const Json& root, const std::string& context)
    -> AtomParameters {
    const auto* atom_json = optional_member(root, "atom");

    if (atom_json == nullptr) {
        return AtomParameters{};
    }

    const auto atom_context = child_context(context, "atom");
    ensure_object(*atom_json, atom_context);

    const auto names = parse_names(*atom_json, atom_context);

    const auto data_context = child_context(atom_context, "data");
    const auto& data_json = required_member(*atom_json, "data", atom_context);

    ensure_array(data_json, data_context);

    std::vector<AtomParameterEntry> entries;
    entries.reserve(data_json.size());

    for (std::size_t entry_index = 0; entry_index < data_json.size(); ++entry_index) {
        const auto entry_context = array_context(data_context, entry_index);
        const auto& entry_json = data_json[entry_index];

        ensure_object(entry_json, entry_context);

        const auto key_context = child_context(entry_context, "key");
        const auto values_context = child_context(entry_context, "values");

        const auto key =
            parse_atom_key(required_member(entry_json, "key", entry_context), key_context);

        auto parameters = make_named_parameters(
            names, required_member(entry_json, "values", entry_context), values_context);

        entries.push_back(AtomParameterEntry{.key = key, .parameters = std::move(parameters)});
    }

    return AtomParameters{std::move(entries)};
}

[[nodiscard]] auto parse_bond_parameters(const Json& root, const std::string& context)
    -> BondParameters {
    const auto* bond_json = optional_member(root, "bond");

    if (bond_json == nullptr) {
        return BondParameters{};
    }

    const auto bond_context = child_context(context, "bond");
    ensure_object(*bond_json, bond_context);

    const auto names = parse_names(*bond_json, bond_context);

    const auto data_context = child_context(bond_context, "data");
    const auto& data_json = required_member(*bond_json, "data", bond_context);

    ensure_array(data_json, data_context);

    std::vector<BondParameterEntry> entries;
    entries.reserve(data_json.size());

    for (std::size_t entry_index = 0; entry_index < data_json.size(); ++entry_index) {
        const auto entry_context = array_context(data_context, entry_index);
        const auto& entry_json = data_json[entry_index];

        ensure_object(entry_json, entry_context);

        const auto key_context = child_context(entry_context, "key");
        const auto values_context = child_context(entry_context, "values");

        const auto key =
            parse_bond_key(required_member(entry_json, "key", entry_context), key_context);

        auto parameters = make_named_parameters(
            names, required_member(entry_json, "values", entry_context), values_context);

        entries.push_back(BondParameterEntry{.key = key, .parameters = std::move(parameters)});
    }

    return BondParameters{std::move(entries)};
}

[[nodiscard]] auto parse_parameter_set(const Json& root, const std::string& context)
    -> ParameterSet {
    ensure_object(root, context);

    return ParameterSet{parse_metadata(root, context), parse_common_parameters(root, context),
                        parse_atom_parameters(root, context), parse_bond_parameters(root, context)};
}

[[nodiscard]] auto load_parameter_set_json_document(std::istream& input, const std::string& context)
    -> ParameterSet {
    try {
        const auto root = Json::parse(input);
        return parse_parameter_set(root, context);
    } catch (const nlohmann::json::exception& error) {
        throw_error(context, error.what());
    }
}

[[nodiscard]] auto json_files_in_directory(const std::filesystem::path& directory)
    -> std::vector<std::filesystem::path> {
    if (!std::filesystem::exists(directory)) {
        throw std::invalid_argument{"parameter directory '" + directory.string() +
                                    "' does not exist"};
    }

    if (!std::filesystem::is_directory(directory)) {
        throw std::invalid_argument{"parameter path '" + directory.string() +
                                    "' is not a directory"};
    }

    std::vector<std::filesystem::path> files;

    for (const auto& entry : std::filesystem::directory_iterator{directory}) {
        if (!entry.is_regular_file()) {
            continue;
        }

        if (entry.path().extension() == ".json") {
            files.push_back(entry.path());
        }
    }

    std::ranges::sort(files);

    return files;
}

auto reject_duplicate_parameter_set_id(const std::set<std::string>& ids,
                                       const ParameterSet& parameter_set,
                                       const std::filesystem::path& source) -> void {
    if (ids.contains(std::string{parameter_set.id()})) {
        throw std::invalid_argument{"duplicate parameter set id '" +
                                    std::string{parameter_set.id()} + "' while loading '" +
                                    source.string() + "'"};
    }
}

auto append_parameter_sets_from_directory_if_present(std::vector<ParameterSet>& parameter_sets,
                                                     std::set<std::string>& ids,
                                                     const std::filesystem::path& directory)
    -> void {
    if (!std::filesystem::exists(directory)) {
        return;
    }

    if (!std::filesystem::is_directory(directory)) {
        throw std::invalid_argument{"default parameter path '" + directory.string() +
                                    "' is not a directory"};
    }

    for (const auto& file : json_files_in_directory(directory)) {
        auto parameter_set = load_parameter_set_json_file(file);
        reject_duplicate_parameter_set_id(ids, parameter_set, file);
        ids.insert(std::string{parameter_set.id()});
        parameter_sets.push_back(std::move(parameter_set));
    }
}

// Group parameter sets by method, each method's sets in the order automatic selection tries them.
auto sort_parameter_sets(std::vector<ParameterSet>& parameter_sets) -> void {
    std::ranges::sort(parameter_sets, [](const ParameterSet& first, const ParameterSet& second) {
        if (first.method_id() != second.method_id()) {
            return first.method_id() < second.method_id();
        }
        return first.priority() != second.priority() ? first.priority() > second.priority()
                                                     : first.id() < second.id();
    });
}

} // namespace

auto load_parameter_set_json(std::istream& input) -> ParameterSet {
    return load_parameter_set_json_document(input, "parameter JSON");
}

auto load_parameter_set_json_file(const std::filesystem::path& path) -> ParameterSet {
    std::ifstream input{path};

    if (!input) {
        throw std::invalid_argument{"failed to open parameter file '" + path.string() + "'"};
    }

    return load_parameter_set_json_document(input, "parameter file '" + path.string() + "'");
}

auto load_parameter_sets_json_directory(const std::filesystem::path& directory)
    -> std::vector<ParameterSet> {
    std::vector<ParameterSet> parameter_sets;
    std::set<std::string> ids;

    for (const auto& file : json_files_in_directory(directory)) {
        auto parameter_set = load_parameter_set_json_file(file);
        reject_duplicate_parameter_set_id(ids, parameter_set, file);
        ids.insert(std::string{parameter_set.id()});
        parameter_sets.push_back(std::move(parameter_set));
    }

    sort_parameter_sets(parameter_sets);
    return parameter_sets;
}

auto load_default_parameter_sets() -> std::vector<ParameterSet> {
    std::vector<ParameterSet> parameter_sets;
    std::set<std::string> ids;

    for (const auto& directory : default_parameter_directories()) {
        append_parameter_sets_from_directory_if_present(parameter_sets, ids, directory);
    }

    if (parameter_sets.empty()) {
        auto searched = std::string{};
        for (const auto& directory : default_parameter_directories()) {
            searched += (searched.empty() ? "'" : ", '") + directory.string() + "'";
        }
        throw std::invalid_argument{"no bundled parameter sets were found in " + searched};
    }

    sort_parameter_sets(parameter_sets);
    return parameter_sets;
}

} // namespace chargefw::parameters
