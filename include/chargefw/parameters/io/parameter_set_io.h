#pragma once

#include <chargefw/parameters/models/parameter_set.h>

#include <filesystem>
#include <istream>
#include <ostream>
#include <string>
#include <vector>

namespace chargefw::parameters {

// Reads one parameter set. Content the reader refuses throws ParseError, which names the source
// when one is given; the file loader gives the path. A file that cannot be opened throws
// std::runtime_error.
[[nodiscard]] auto load_parameter_set_json(std::istream& input, std::string source_name = {})
    -> ParameterSet;

[[nodiscard]] auto load_parameter_set_json_file(const std::filesystem::path& path) -> ParameterSet;

// Writes `set` in the parameter-set JSON format, which load_parameter_set_json() reads back as an
// equal set. The format lists parameter names once per section, so every atom entry, and every bond
// entry, must have the same names in the same order. Keys must be ones the reader accepts: a known
// element or the wildcard, and type "*" for a plain classifier. Otherwise std::invalid_argument,
// and nothing is written. A section without entries is omitted.
auto save_parameter_set_json(const ParameterSet& set, std::ostream& output) -> void;

// The sets of the .json files in `directory`. A directory that does not exist and two files with
// the same parameter-set ID throw std::invalid_argument.
[[nodiscard]] auto load_parameter_sets_json_directory(const std::filesystem::path& directory)
    -> std::vector<ParameterSet>;

// The bundled sets. An installation without them, or with two sets of the same ID, throws
// std::runtime_error.
[[nodiscard]] auto load_default_parameter_sets() -> std::vector<ParameterSet>;

} // namespace chargefw::parameters