#pragma once

#include <cstddef>

namespace chargefw::platform {

// Physical processor cores this process may run on, or 0 when the platform does not report them.
// Simultaneous threads of one core count once.
[[nodiscard]] auto available_physical_cores() -> std::size_t;

} // namespace chargefw::platform
