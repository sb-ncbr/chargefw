#pragma once

#include <chargefw/methods/method.h>

namespace chargefw::methods::builtin {

class FormalMethod final : public Method {
  public:
    [[nodiscard]] auto metadata() const noexcept -> const MethodMetadata& override {
        static constexpr MethodMetadata metadata{
            .id = "formal",
            .name = "Formal",
            .full_name = "Formal atomic charges",
            .publication = std::nullopt,
            .notes = "Copies input atomic formal charges unchanged; "
                     "not an empirical charge model.",
            .priority = 10,
            .recommended = false};

        return metadata;
    }

    [[nodiscard]] auto requirements() const noexcept -> MethodRequirements override {
        auto requirements = MethodRequirements{};
        requirements.resources.time = ComplexityTerm::atoms;
        requirements.resources.memory = ComplexityTerm::constant;
        return requirements;
    }

    [[nodiscard]] auto option_schema() const noexcept
        -> std::span<const MethodOptionSpec> override {
        return {};
    }

    [[nodiscard]] auto calculate(const CalculationInput& input) const
        -> charges::AtomicCharges override;

    [[nodiscard]] auto estimate(const TargetSize& target, const MethodOptions& options) const
        -> ResourceEstimate override;
};

} // namespace chargefw::methods::builtin
