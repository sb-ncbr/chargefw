#pragma once

#include <cstdint>

namespace chargefw::io::native {

enum class BondFormat : std::uint8_t {
    mol,
    mol2,
};

} // namespace chargefw::io::native
