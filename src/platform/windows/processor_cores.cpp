#include "platform/processor_cores.h"

#include <vector>

// Included last so that its macros cannot affect other headers.
#define WIN32_LEAN_AND_MEAN
#include <windows.h>

namespace chargefw::platform {

// Counts the cores of the machine; a process affinity mask is not applied.
auto available_physical_cores() -> std::size_t {
    DWORD length = 0;
    GetLogicalProcessorInformationEx(RelationProcessorCore, nullptr, &length);
    if (length == 0) {
        return 0;
    }

    auto buffer = std::vector<char>(length);
    if (GetLogicalProcessorInformationEx(
            RelationProcessorCore,
            reinterpret_cast<SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX*>(buffer.data()),
            &length) == 0) {
        return 0;
    }

    // One variable-size record per core.
    std::size_t cores = 0;
    for (DWORD offset = 0; offset < length; ++cores) {
        offset +=
            reinterpret_cast<const SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX*>(buffer.data() + offset)
                ->Size;
    }
    return cores;
}

} // namespace chargefw::platform
