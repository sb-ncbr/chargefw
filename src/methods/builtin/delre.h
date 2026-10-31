#pragma once

#include <chargefw/methods/method.h>

namespace chargefw::methods::builtin {

class DelReMethod final : public Method {
  public:
    [[nodiscard]] auto metadata() const noexcept -> const MethodMetadata& override {
        static constexpr MethodMetadata metadata{.id = "delre",
                                                 .name = "DelRe",
                                                 .full_name = "Method of Del Re",
                                                 .publication = "10.1039/JR9580004031",
                                                 .notes =
                                                     "Every connected component must be formally "
                                                     "neutral.",
                                                 .priority = 130,
                                                 .recommended = false};

        return metadata;
    }

    [[nodiscard]] auto requirements() const -> MethodRequirements override {
        auto requirements = MethodRequirements{};
        requirements.bonds = true;
        requirements.atom_parameters = {"delta"};
        requirements.bond_parameters = {"eps", "gammaA", "gammaB"};
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
