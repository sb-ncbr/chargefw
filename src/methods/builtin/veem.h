#pragma once

#include <chargefw/methods/method.h>

namespace chargefw::methods::builtin {

class VEEMMethod final : public Method {
  public:
    [[nodiscard]] auto metadata() const noexcept -> const MethodMetadata& override {
        static constexpr MethodMetadata metadata{
            .id = "veem",
            .name = "VEEM",
            .full_name = "Valence Electrons Equalization Method",
            .publication = "10.1088/1674-0068/24/01/31-39",
            .notes = "Supports only net-neutral molecules; molecular "
                     "connectivity is not used.",
            .priority = 20,
            .recommended = false};

        return metadata;
    }

    [[nodiscard]] auto requirements() const -> MethodRequirements override {
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

  protected:
    auto add_method_specific_prerequisite_issues(const MethodPrerequisiteInput& input,
                                                 PrerequisiteResult& result) const -> void override;
};

} // namespace chargefw::methods::builtin
