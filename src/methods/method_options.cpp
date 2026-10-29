#include <chargefw/methods/method_options.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <utility>
#include <vector>

namespace chargefw::methods {

auto to_string(const MethodOptionType value) -> std::string_view {
    switch (value) {
    case MethodOptionType::boolean:
        return "boolean";
    case MethodOptionType::integer:
        return "integer";
    case MethodOptionType::floating_point:
        return "floating_point";
    case MethodOptionType::string:
        return "string";
    }

    throw std::invalid_argument{"unknown method option type"};
}

namespace {

[[nodiscard]] auto value_type(const MethodOptionValue& value) -> MethodOptionType {
    if (std::holds_alternative<bool>(value)) {
        return MethodOptionType::boolean;
    }

    if (std::holds_alternative<int>(value)) {
        return MethodOptionType::integer;
    }

    if (std::holds_alternative<double>(value)) {
        return MethodOptionType::floating_point;
    }

    return MethodOptionType::string;
}

[[nodiscard]] auto value_matches_type(const MethodOptionValue& value,
                                      const MethodOptionType expected_type) -> bool {
    return value_type(value) == expected_type;
}

[[nodiscard]] auto choices_contain(const std::vector<MethodOptionValue>& choices,
                                   const MethodOptionValue& value) -> bool {
    return std::ranges::find(choices, value) != choices.end();
}

[[nodiscard]] auto numeric_value(const MethodOptionValue& value) -> double {
    if (const auto* integer = std::get_if<int>(&value); integer != nullptr) {
        return static_cast<double>(*integer);
    }

    return std::get<double>(value);
}

auto validate_bounds(const MethodOptionSpec& option, const MethodOptionValue& value) -> void {
    // NaN fails every comparison and infinity passes a lone lower bound, so check finiteness first.
    if (const auto* number = std::get_if<double>(&value);
        number != nullptr && !std::isfinite(*number)) {
        throw std::invalid_argument{"method option '" + std::string{option.id} +
                                    "' must be a finite number"};
    }
    if (!option.minimum.has_value() && !option.maximum.has_value()) {
        return;
    }

    const auto number = numeric_value(value);
    if (option.minimum.has_value() &&
        (number < numeric_value(*option.minimum) ||
         (!option.minimum_inclusive && number == numeric_value(*option.minimum)))) {
        throw std::invalid_argument{"method option '" + std::string{option.id} + "' must be " +
                                    (option.minimum_inclusive ? "at least " : "greater than ") +
                                    std::to_string(numeric_value(*option.minimum))};
    }
    if (option.maximum.has_value() &&
        (number > numeric_value(*option.maximum) ||
         (!option.maximum_inclusive && number == numeric_value(*option.maximum)))) {
        throw std::invalid_argument{"method option '" + std::string{option.id} + "' must be " +
                                    (option.maximum_inclusive ? "at most " : "less than ") +
                                    std::to_string(numeric_value(*option.maximum))};
    }
}

[[nodiscard]] auto find_spec(std::span<const MethodOptionSpec> schema, const std::string_view id)
    -> const MethodOptionSpec* {
    const auto iter = std::ranges::find_if(
        schema, [id](const MethodOptionSpec& option) -> bool { return option.id == id; });

    if (iter == schema.end()) {
        return nullptr;
    }

    return std::addressof(*iter);
}

} // namespace

MethodOptions::MethodOptions(std::unordered_map<std::string, MethodOptionValue> values)
    : values_{std::move(values)} {}

auto MethodOptions::contains(const std::string_view id) const -> bool {
    return values_.contains(std::string{id});
}

auto MethodOptions::values() const noexcept
    -> const std::unordered_map<std::string, MethodOptionValue>& {
    return values_;
}

auto MethodOptions::set(std::string id, MethodOptionValue value) -> void {
    values_.insert_or_assign(std::move(id), std::move(value));
}

auto validate_method_option_schema(std::span<const MethodOptionSpec> schema) -> void {
    std::vector<std::string_view> ids;
    ids.reserve(schema.size());

    for (const auto& option : schema) {
        if (option.id.empty()) {
            throw std::invalid_argument{"method option schema contains an empty option id"};
        }

        if (std::ranges::find(ids, option.id) != ids.end()) {
            throw std::invalid_argument{"duplicate method option id '" + std::string{option.id} +
                                        "'"};
        }

        ids.push_back(option.id);

        if (!value_matches_type(option.default_value, option.type)) {
            throw std::invalid_argument{
                "default value for method option '" + std::string{option.id} + "' has type '" +
                std::string{to_string(value_type(option.default_value))} + "', expected '" +
                std::string{to_string(option.type)} + "'"};
        }

        if ((option.minimum.has_value() || option.maximum.has_value()) &&
            option.type != MethodOptionType::integer &&
            option.type != MethodOptionType::floating_point) {
            throw std::invalid_argument{"method option '" + std::string{option.id} +
                                        "' has numeric bounds but is not numeric"};
        }
        if ((option.minimum.has_value() && !value_matches_type(*option.minimum, option.type)) ||
            (option.maximum.has_value() && !value_matches_type(*option.maximum, option.type))) {
            throw std::invalid_argument{"method option '" + std::string{option.id} +
                                        "' has a bound with an unexpected type"};
        }
        if ((option.minimum.has_value() && !std::isfinite(numeric_value(*option.minimum))) ||
            (option.maximum.has_value() && !std::isfinite(numeric_value(*option.maximum)))) {
            throw std::invalid_argument{"method option '" + std::string{option.id} +
                                        "' has a non-finite bound"};
        }
        if (option.minimum.has_value() && option.maximum.has_value() &&
            (numeric_value(*option.minimum) > numeric_value(*option.maximum) ||
             (numeric_value(*option.minimum) == numeric_value(*option.maximum) &&
              (!option.minimum_inclusive || !option.maximum_inclusive)))) {
            throw std::invalid_argument{"method option '" + std::string{option.id} +
                                        "' has an empty numeric range"};
        }
        validate_bounds(option, option.default_value);

        for (const auto& choice : option.choices) {
            if (!value_matches_type(choice, option.type)) {
                throw std::invalid_argument{
                    "choice value for method option '" + std::string{option.id} + "' has type '" +
                    std::string{to_string(value_type(choice))} + "', expected '" +
                    std::string{to_string(option.type)} + "'"};
            }
            validate_bounds(option, choice);
        }

        if (!option.choices.empty() && !choices_contain(option.choices, option.default_value)) {
            throw std::invalid_argument{"default value for method option '" +
                                        std::string{option.id} +
                                        "' is not one of the allowed choices"};
        }
    }
}

auto make_default_options(std::span<const MethodOptionSpec> schema) -> MethodOptions {
    validate_method_option_schema(schema);

    MethodOptions options;

    for (const auto& option : schema) {
        options.set(std::string{option.id}, option.default_value);
    }

    return options;
}

auto validate_method_options(std::span<const MethodOptionSpec> schema, const MethodOptions& options)
    -> void {
    validate_method_option_schema(schema);

    for (const auto& [id, value] : options.values()) {
        const auto* spec = find_spec(schema, id);

        if (spec == nullptr) {
            throw std::invalid_argument{"unknown method option '" + id + "'"};
        }

        if (!value_matches_type(value, spec->type)) {
            throw std::invalid_argument{"method option '" + id + "' has type '" +
                                        std::string{to_string(value_type(value))} +
                                        "', expected '" + std::string{to_string(spec->type)} + "'"};
        }

        if (!spec->choices.empty() && !choices_contain(spec->choices, value)) {
            throw std::invalid_argument{"method option '" + id +
                                        "' is not one of the allowed choices"};
        }
        validate_bounds(*spec, value);
    }

    for (const auto& option : schema) {
        if (!options.contains(option.id)) {
            throw std::invalid_argument{"missing method option '" + std::string{option.id} + "'"};
        }
    }
}

} // namespace chargefw::methods
