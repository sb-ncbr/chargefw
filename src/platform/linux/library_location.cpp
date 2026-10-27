#include "platform/library_location.h"

#include <dlfcn.h>
#include <link.h>
#include <stdexcept>

namespace chargefw::platform {

auto loaded_library_path() -> std::filesystem::path {
    Dl_info information{};
    link_map* map = nullptr;
    if (dladdr1(reinterpret_cast<const void*>(&loaded_library_path), &information,
                reinterpret_cast<void**>(&map), RTLD_DL_LINKMAP) == 0 ||
        map == nullptr) {
        throw std::runtime_error{"could not determine the ChargeFW library location"};
    }

    // For the main program, glibc reports its argv[0], which need not be a path to the file; its
    // link-map entry has an empty name instead.
    if (map->l_name[0] == '\0') {
        return std::filesystem::read_symlink("/proc/self/exe");
    }

    return std::filesystem::absolute(map->l_name);
}

} // namespace chargefw::platform
