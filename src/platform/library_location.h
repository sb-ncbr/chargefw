#pragma once

#include <filesystem>

namespace chargefw::platform {

// The file of the loaded library, or of the executable or extension module that ChargeFW is linked
// into statically.
[[nodiscard]] auto loaded_library_path() -> std::filesystem::path;

} // namespace chargefw::platform
