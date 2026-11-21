#include "platform/processor_cores.h"

#include <fstream>
#include <set>
#include <string>

#include <sched.h>

namespace chargefw::platform {

auto available_physical_cores() -> std::size_t {
    cpu_set_t allowed;
    CPU_ZERO(&allowed);
    if (sched_getaffinity(0, sizeof(allowed), &allowed) != 0) {
        return 0;
    }

    // The processors of one core list the same siblings, such as "0,8" or "0-1".
    auto cores = std::set<std::string>{};
    for (std::size_t processor = 0; processor < CPU_SETSIZE; ++processor) {
        if (!CPU_ISSET(processor, &allowed)) {
            continue;
        }
        auto file = std::ifstream{"/sys/devices/system/cpu/cpu" + std::to_string(processor) +
                                  "/topology/thread_siblings_list"};
        auto siblings = std::string{};
        if (!std::getline(file, siblings)) {
            return 0;
        }
        cores.insert(std::move(siblings));
    }
    return cores.size();
}

} // namespace chargefw::platform
