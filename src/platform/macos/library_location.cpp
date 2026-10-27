#include "platform/library_location.h"

#include <dlfcn.h>
#include <stdexcept>

namespace chargefw::platform {

auto loaded_library_path() -> std::filesystem::path {
    Dl_info information{};
    if (dladdr(reinterpret_cast<const void*>(&loaded_library_path), &information) == 0 ||
        information.dli_fname == nullptr) {
        throw std::runtime_error{"could not determine the ChargeFW library location"};
    }

    return std::filesystem::absolute(information.dli_fname);
}

} // namespace chargefw::platform
