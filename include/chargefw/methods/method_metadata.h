#pragma once

#include <optional>
#include <string_view>

namespace chargefw::methods {
struct MethodMetadata {
    std::string_view id;
    std::string_view name;
    std::string_view full_name;
    std::optional<std::string_view> publication;
    std::string_view notes{};
    int priority = 0;
    // Whether the method's charges are recommended for general use. Automatic selection considers
    // only recommended methods unless the request includes all methods.
    bool recommended = true;
};

} // namespace chargefw::methods
