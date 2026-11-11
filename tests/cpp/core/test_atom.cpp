#include <chargefw/core/atom.h>

#include <stdexcept>

#include <snitch/snitch.hpp>

namespace core = chargefw::core;

TEST_CASE("atom rejects invalid atomic numbers", "[core][atom]") {
    CHECK(core::Atom{100}.atomic_number() == 100);
    CHECK_THROWS_AS(core::Atom{0}, std::invalid_argument);
    CHECK_THROWS_AS(core::Atom{101}, std::invalid_argument);
    CHECK_THROWS_AS(core::Atom{118}, std::invalid_argument);
    CHECK_THROWS_AS(core::Atom{119}, std::invalid_argument);
}
