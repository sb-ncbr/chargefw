#pragma once

#include <chargefw/parameters/models/parameter_set.h>

#include <filesystem>
#include <istream>
#include <ostream>
#include <vector>

namespace chargefw::parameters {

[[nodiscard]] auto load_parameter_set_json(std::istream& input) -> ParameterSet;

[[nodiscard]] auto load_parameter_set_json_file(const std::filesystem::path& path) -> ParameterSet;

// Writes `set` in the parameter-set JSON format, which load_parameter_set_json() reads back as an
// equal set. The format lists parameter names once per section, so every atom entry, and every bond
// entry, must have the same names in the same order. Keys must be ones the reader accepts: a known
// element or the wildcard, and type "*" for a plain classifier. Otherwise std::invalid_argument,
// and nothing is written. A section without entries is omitted.
auto save_parameter_set_json(const ParameterSet& set, std::ostream& output) -> void;

[[nodiscard]] auto load_parameter_sets_json_directory(const std::filesystem::path& directory)
    -> std::vector<ParameterSet>;

[[nodiscard]] auto load_default_parameter_sets() -> std::vector<ParameterSet>;

} // namespace chargefw::parameters