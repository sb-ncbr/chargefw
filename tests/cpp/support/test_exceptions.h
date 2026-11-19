#pragma once

#include <optional>
#include <utility>

namespace chargefw::test {

// The exception of type Error that `call` throws, so that a test can check its fields; nothing
// when `call` throws none. Other exceptions pass.
template <typename Error, typename Call>
[[nodiscard]] auto thrown(Call&& call) -> std::optional<Error> {
    try {
        std::forward<Call>(call)();
    } catch (const Error& error) {
        return error;
    }
    return std::nullopt;
}

} // namespace chargefw::test
