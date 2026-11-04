#pragma once

#include <chargefw/methods/method.h>

namespace chargefw::methods::builtin {

class EQeqCMethod final : public Method {
  public:
    [[nodiscard]] auto metadata() const noexcept -> const MethodMetadata& override {
        static constexpr MethodMetadata metadata{
            .id = "eqeqc",
            .name = "EQeq+C",
            .full_name = "Bond-Order-Corrected Extended Charge Equilibration Method",
            .publication = "10.1021/acs.jctc.5b00037",
            .notes = "The correction depends on the distances between atoms; the bonds of the "
                     "input are not used.",
            .priority = 140};

        return metadata;
    }

    [[nodiscard]] auto requirements() const -> MethodRequirements override {
        auto requirements = MethodRequirements{};
        requirements.coordinates = true;
        requirements.common_parameters = {"alpha"};
        requirements.atom_parameters = {"Dz"};
        requirements.resources.time = ComplexityTerm::atoms_cubed;
        requirements.resources.memory = ComplexityTerm::atoms_squared;
        requirements.resources.supports_cutoff = true;
        requirements.resources.supports_cover = true;
        requirements.resources.reduced_charge_policy = ReducedChargePolicy::uniform_target_global;
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
