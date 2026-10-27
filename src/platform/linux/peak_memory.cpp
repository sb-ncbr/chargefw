#include "platform/peak_memory.h"

#include <sys/resource.h>

namespace chargefw::platform {

auto peak_resident_memory_bytes() -> std::uint64_t {
    rusage usage{};
    if (getrusage(RUSAGE_SELF, &usage) != 0) {
        return 0;
    }

    // Linux reports ru_maxrss in kilobytes.
    return static_cast<std::uint64_t>(usage.ru_maxrss) * 1024;
}

} // namespace chargefw::platform
