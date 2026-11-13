#include <chargefw/calculation/execution_policy.h>

#include <cmath>
#include <stdexcept>
#include <string>

namespace chargefw::calculation {
namespace {

auto validate_reduced_radius(const std::optional<double> radius) -> void {
    if (!radius.has_value()) {
        throw std::invalid_argument{"cutoff and cover execution require a radius"};
    }

    if (!std::isfinite(*radius)) {
        throw std::invalid_argument{"cutoff and cover radius must be finite"};
    }

    if (*radius < minimum_reduced_radius) {
        throw std::invalid_argument{"cutoff and cover radius must be at least " +
                                    std::to_string(minimum_reduced_radius) + " angstrom"};
    }
}

auto validate_policy(const ExecutionMode mode, const std::optional<double> radius) -> void {
    switch (mode) {
    case ExecutionMode::full:
        if (radius.has_value()) {
            throw std::invalid_argument{"full execution does not accept a radius"};
        }
        return;
    case ExecutionMode::cutoff:
    case ExecutionMode::cover:
        validate_reduced_radius(radius);
        return;
    }

    throw std::invalid_argument{"unknown execution mode"};
}

auto validate_selection(const ExecutionSelectionMode kind, const std::optional<double> radius)
    -> void {
    switch (kind) {
    case ExecutionSelectionMode::automatic:
        if (radius.has_value()) {
            validate_reduced_radius(radius);
        }
        return;
    case ExecutionSelectionMode::full:
        validate_policy(ExecutionMode::full, radius);
        return;
    case ExecutionSelectionMode::cutoff:
    case ExecutionSelectionMode::cover:
        validate_reduced_radius(radius);
        return;
    }

    throw std::invalid_argument{"unknown execution selection"};
}

} // namespace

auto execution_selection_mode_from_string(const std::string_view value) -> ExecutionSelectionMode {
    if (value == "auto") {
        return ExecutionSelectionMode::automatic;
    }
    if (value == "full") {
        return ExecutionSelectionMode::full;
    }
    if (value == "cutoff") {
        return ExecutionSelectionMode::cutoff;
    }
    if (value == "cover") {
        return ExecutionSelectionMode::cover;
    }
    throw std::invalid_argument{"unknown execution selection: " + std::string{value}};
}

auto to_string(const ExecutionSelectionMode value) -> std::string_view {
    switch (value) {
    case ExecutionSelectionMode::automatic:
        return "auto";
    case ExecutionSelectionMode::full:
        return "full";
    case ExecutionSelectionMode::cutoff:
        return "cutoff";
    case ExecutionSelectionMode::cover:
        return "cover";
    }
    throw std::invalid_argument{"unknown execution selection"};
}

auto to_string(const ExecutionMode value) -> std::string_view {
    switch (value) {
    case ExecutionMode::full:
        return "full";
    case ExecutionMode::cutoff:
        return "cutoff";
    case ExecutionMode::cover:
        return "cover";
    }
    throw std::invalid_argument{"unknown execution mode"};
}

ExecutionPolicy::ExecutionPolicy(const ExecutionMode mode, const std::optional<double> radius)
    : mode_{mode}, radius_{radius} {
    validate_policy(mode_, radius_);
}

auto ExecutionPolicy::mode() const noexcept -> ExecutionMode {
    return mode_;
}

auto ExecutionPolicy::radius() const noexcept -> std::optional<double> {
    return radius_;
}

ExecutionSelection::ExecutionSelection(const ExecutionSelectionMode mode,
                                       const std::optional<double> radius)
    : mode_{mode}, radius_{radius} {
    validate_selection(mode_, radius_);
}

auto ExecutionSelection::mode() const noexcept -> ExecutionSelectionMode {
    return mode_;
}

auto ExecutionSelection::radius() const noexcept -> std::optional<double> {
    return radius_;
}

} // namespace chargefw::calculation
