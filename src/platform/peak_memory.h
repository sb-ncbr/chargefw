#pragma once

#include <cstdint>

namespace chargefw::platform {

// Peak resident memory of this process in bytes, or 0 when the platform does not report it.
[[nodiscard]] auto peak_resident_memory_bytes() -> std::uint64_t;

} // namespace chargefw::platform
