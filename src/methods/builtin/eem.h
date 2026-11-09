#pragma once

#include <chargefw/methods/method.h>

namespace chargefw::methods::builtin {

class EEMMethod final : public Method {
  public:
    [[nodiscard]] auto metadata() const noexcept -> const MethodMetadata& override {
        static constexpr MethodMetadata metadata{.id = "eem",
                                                 .name = "EEM",
                                                 .full_name =
                                                     "Electronegativity Equalization Method",
                                                 .publication = "10.1021/ja00275a013",
                                                 .priority = 200};

        return metadata;
    }

    [[nodiscard]] auto requirements() const -> MethodRequirements override {
        auto requirements = MethodRequirements{};
        requirements.coordinates = true;
        requirements.common_parameters = {"kappa"};
        requirements.atom_parameters = {"A", "B"};
        requirements.resources.time = ComplexityTerm::atoms_cubed;
        requirements.resources.memory = ComplexityTerm::atoms_squared;
        requirements.resources.supports_cutoff = true;
        requirements.resources.supports_cover = true;
        requirements.resources.reduced_charge_policy = ReducedChargePolicy::uniform_target_global;
        requirements.supports_fixed_point_sources = true;
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
