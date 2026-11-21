#include "platform/processor_cores.h"

#include <sys/sysctl.h>
#include <sys/types.h>

namespace chargefw::platform {

auto available_physical_cores() -> std::size_t {
    int cores = 0;
    auto size = sizeof(cores);
    if (sysctlbyname("hw.physicalcpu", &cores, &size, nullptr, 0) != 0 || cores < 0) {
        return 0;
    }
    return static_cast<std::size_t>(cores);
}

} // namespace chargefw::platform
