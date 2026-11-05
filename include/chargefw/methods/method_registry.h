#pragma once

#include <chargefw/methods/method.h>

#include <memory>
#include <span>
#include <string_view>
#include <vector>

namespace chargefw::methods {

// The read-only catalog of built-in methods, available through method_registry().
class MethodRegistry {
  public:
    [[nodiscard]] auto find(std::string_view id) const noexcept -> const Method*;

    [[nodiscard]] auto methods() const noexcept -> std::span<const std::unique_ptr<Method>>;

  private:
    explicit MethodRegistry(std::vector<std::unique_ptr<Method>> methods);

    friend auto method_registry() -> const MethodRegistry&;

    std::vector<std::unique_ptr<Method>> methods_;
};

[[nodiscard]] auto method_registry() -> const MethodRegistry&;

} // namespace chargefw::methods
