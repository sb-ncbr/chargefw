#pragma once

#include <chargefw/methods/method.h>

namespace chargefw::methods::builtin {

class SMPQEqMethod final : public Method {
  public:
    [[nodiscard]] auto metadata() const noexcept -> const MethodMetadata& override {
        static constexpr MethodMetadata metadata{
            .id = "smpqeq",
            .name = "SMP/QEq",
            .full_name = "Self-Consistent Charge Equilibration Method",
            .publication = "10.1021/jp8063273",
            .notes = "Iterates charge-dependent hardness to convergence.",
            .priority = 160,
            .recommended = false};

        return metadata;
    }

    [[nodiscard]] auto requirements() const -> MethodRequirements override {
        auto requirements = MethodRequirements{};
        requirements.coordinates = true;
        requirements.atom_parameters = {"first", "second", "third", "fourth"};
        requirements.resources.time = ComplexityTerm::atoms_cubed;
        requirements.resources.memory = ComplexityTerm::atoms_squared;
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
