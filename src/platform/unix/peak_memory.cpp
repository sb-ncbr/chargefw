#include "platform/peak_memory.h"

#include <sys/resource.h>

namespace chargefw::platform {

auto peak_resident_memory_bytes() -> std::uint64_t {
    rusage usage{};
    if (getrusage(RUSAGE_SELF, &usage) != 0) {
        return 0;
    }

    // POSIX leaves the unit of ru_maxrss unspecified: Linux reports kilobytes, macOS bytes.
#if defined(__APPLE__)
    constexpr std::uint64_t unit = 1;
#else
    constexpr std::uint64_t unit = 1024;
#endif
    return static_cast<std::uint64_t>(usage.ru_maxrss) * unit;
}

} // namespace chargefw::platform
