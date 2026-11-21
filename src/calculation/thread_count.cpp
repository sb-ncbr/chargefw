#include "calculation/parallel_for.h"
#include "platform/processor_cores.h"

#include <algorithm>
#include <cstddef>

#include <oneapi/tbb/global_control.h>

namespace chargefw::calculation::detail {

auto thread_count(const std::size_t max_threads) -> std::size_t {
    const auto allowed = oneapi::tbb::global_control::active_value(
        oneapi::tbb::global_control::max_allowed_parallelism);
    if (max_threads != 0) {
        return std::min(max_threads, allowed);
    }
    // Read once: the processors of the machine do not change while the process runs.
    static const auto cores = platform::available_physical_cores();
    return cores == 0 ? allowed : std::min(cores, allowed);
}

} // namespace chargefw::calculation::detail
