#include <chargefw/parameters/io/parameter_data_paths.h>

#include "platform/library_location.h"

#include <chargefw/config.h>

#include <filesystem>
#include <vector>

namespace chargefw::parameters {

auto bundled_parameter_directory() -> std::filesystem::path {
    return (platform::loaded_library_path().parent_path() /
            CHARGEFW_PARAMETER_DIRECTORY_FROM_LIBRARY)
        .lexically_normal();
}

auto default_parameter_directories() -> std::vector<std::filesystem::path> {
    return {bundled_parameter_directory()};
}

} // namespace chargefw::parameters
