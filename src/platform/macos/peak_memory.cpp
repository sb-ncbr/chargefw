#include "platform/peak_memory.h"

#include <sys/resource.h>

namespace chargefw::platform {

auto peak_resident_memory_bytes() -> std::uint64_t {
    rusage usage{};
    if (getrusage(RUSAGE_SELF, &usage) != 0) {
        return 0;
    }

    // macOS reports ru_maxrss in bytes.
    return static_cast<std::uint64_t>(usage.ru_maxrss);
}

} // namespace chargefw::platform
