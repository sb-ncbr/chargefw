#include "platform/peak_memory.h"

// Included last so that their macros cannot affect other headers.
#define WIN32_LEAN_AND_MEAN
#include <windows.h>

#include <psapi.h>

namespace chargefw::platform {

auto peak_resident_memory_bytes() -> std::uint64_t {
    PROCESS_MEMORY_COUNTERS counters{};
    if (GetProcessMemoryInfo(GetCurrentProcess(), &counters, sizeof(counters)) == 0) {
        return 0;
    }

    return counters.PeakWorkingSetSize;
}

} // namespace chargefw::platform
