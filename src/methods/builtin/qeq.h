#pragma once

#include <chargefw/methods/method.h>

#include <array>

namespace chargefw::methods::builtin {

class QEqMethod final : public Method {
  public:
    [[nodiscard]] auto metadata() const noexcept -> const MethodMetadata& override {
        static constexpr MethodMetadata metadata{
            .id = "qeq",
            .name = "QEq",
            .full_name = "Charge Equilibration",
            .publication = "10.1021/j100161a070",
            .notes = "The hardness of hydrogen does not depend on its charge, so the equations are "
                     "solved once; overlap_term selects an empirical Coulomb term.",
            .priority = 170};

        return metadata;
    }

    [[nodiscard]] auto requirements() const -> MethodRequirements override {
        auto requirements = MethodRequirements{};
        requirements.coordinates = true;
        requirements.atom_parameters = {"electronegativity", "hardness"};
        requirements.resources.time = ComplexityTerm::atoms_cubed;
        requirements.resources.memory = ComplexityTerm::atoms_squared;
        requirements.resources.supports_cutoff = true;
        requirements.resources.supports_cover = true;
        requirements.resources.reduced_charge_policy = ReducedChargePolicy::uniform_target_global;
        return requirements;
    }

    [[nodiscard]] auto option_schema() const noexcept
        -> std::span<const MethodOptionSpec> override {
        static const std::array option_schema{MethodOptionSpec{
            .id = "overlap_term",
            .description = "QEq empirical two-center Coulomb term",
            .type = MethodOptionType::string,
            .default_value = std::string{"Nishimoto-Mataga"},
            .choices = {std::string{"Nishimoto-Mataga"}, std::string{"Nishimoto-Mataga-Weiss"},
                        std::string{"Ohno"}, std::string{"Ohno-Klopman"},
                        std::string{"DasGupta-Huzinaga"}}}};

        return {option_schema.data(), option_schema.size()};
    }

    [[nodiscard]] auto calculate(const CalculationInput& input) const
        -> charges::AtomicCharges override;

    [[nodiscard]] auto estimate(const TargetSize& target, const MethodOptions& options) const
        -> ResourceEstimate override;
};

} // namespace chargefw::methods::builtin
