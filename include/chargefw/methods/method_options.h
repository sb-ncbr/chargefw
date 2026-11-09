#pragma once

#include <cstdint>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <variant>
#include <vector>

namespace chargefw::methods {

enum class MethodOptionType : std::uint8_t {
    boolean,
    integer,
    floating_point,
    string,
};

[[nodiscard]] auto to_string(MethodOptionType value) -> std::string_view;

using MethodOptionValue = std::variant<bool, int, double, std::string>;

struct MethodOptionSpec {
    std::string_view id;
    std::string_view description;
    MethodOptionType type;
    MethodOptionValue default_value;
    std::vector<MethodOptionValue> choices;
    std::optional<MethodOptionValue> minimum = std::nullopt;
    bool minimum_inclusive = true;
    std::optional<MethodOptionValue> maximum = std::nullopt;
    bool maximum_inclusive = true;
};

class MethodOptions {
  public:
    MethodOptions() = default;

    explicit MethodOptions(std::unordered_map<std::string, MethodOptionValue> values);

    [[nodiscard]] auto contains(std::string_view id) const -> bool;

    [[nodiscard]] auto values() const noexcept
        -> const std::unordered_map<std::string, MethodOptionValue>&;

    template <typename T> [[nodiscard]] auto get(std::string_view id) const -> const T&;

    auto set(std::string id, MethodOptionValue value) -> void;

  private:
    std::unordered_map<std::string, MethodOptionValue> values_;
};

template <typename T> auto MethodOptions::get(const std::string_view id) const -> const T& {
    const auto iter = values_.find(std::string{id});

    if (iter == values_.end()) {
        throw std::out_of_range{"method option '" + std::string{id} + "' was not found"};
    }

    const auto* value = std::get_if<T>(&iter->second);

    if (value == nullptr) {
        throw std::invalid_argument{"method option '" + std::string{id} + "' has unexpected type"};
    }

    return *value;
}

auto validate_method_option_schema(std::span<const MethodOptionSpec> schema) -> void;

[[nodiscard]] auto make_default_options(std::span<const MethodOptionSpec> schema) -> MethodOptions;

auto validate_method_options(std::span<const MethodOptionSpec> schema, const MethodOptions& options)
    -> void;

} // namespace chargefw::methods
