#include "platform/library_location.h"

#include <stdexcept>
#include <string>

// Included last so that its macros cannot affect other headers.
#define WIN32_LEAN_AND_MEAN
#include <windows.h>

namespace chargefw::platform {

auto loaded_library_path() -> std::filesystem::path {
    HMODULE module = nullptr;
    if (GetModuleHandleExW(GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS |
                               GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT,
                           reinterpret_cast<LPCWSTR>(&loaded_library_path), &module) == 0) {
        throw std::runtime_error{"could not determine the ChargeFW library location"};
    }

    // GetModuleFileNameW truncates silently and returns the buffer size when the path does not fit.
    std::wstring path(MAX_PATH, L'\0');
    while (true) {
        const auto length =
            GetModuleFileNameW(module, path.data(), static_cast<DWORD>(path.size()));
        if (length == 0) {
            throw std::runtime_error{"could not determine the ChargeFW library location"};
        }
        if (length < path.size()) {
            path.resize(length);
            return path;
        }
        path.resize(path.size() * 2);
    }
}

} // namespace chargefw::platform
