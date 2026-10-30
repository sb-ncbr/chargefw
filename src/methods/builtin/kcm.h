#pragma once

#include <chargefw/methods/method.h>

namespace chargefw::methods::builtin {

class KCMMethod final : public Method {
  public:
    [[nodiscard]] auto metadata() const noexcept -> const MethodMetadata& override {
        static constexpr MethodMetadata metadata{.id = "kcm",
                                                 .name = "KCM",
                                                 .full_name = "Kirchhoff Charge Model",
                                                 .publication = "10.1002/jcc.20892",
                                                 .notes =
                                                     "Every connected component must be formally "
                                                     "neutral.",
                                                 .priority = 60};

        return metadata;
    }

    [[nodiscard]] auto requirements() const -> MethodRequirements override {
        auto requirements = MethodRequirements{};
        requirements.bonds = true;
        requirements.atom_parameters = {"electronegativity", "hardness"};
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

  protected:
    auto add_method_specific_prerequisite_issues(const MethodPrerequisiteInput& input,
                                                 PrerequisiteResult& result) const -> void override;
};

} // namespace chargefw::methods::builtin
