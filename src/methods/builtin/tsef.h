#pragma once

#include <chargefw/methods/method.h>

namespace chargefw::methods::builtin {

class TSEFMethod final : public Method {
  public:
    [[nodiscard]] auto metadata() const noexcept -> const MethodMetadata& override {
        static constexpr MethodMetadata metadata{
            .id = "tsef",
            .name = "TSEF",
            .full_name = "Topologically Symmetrical Energy Function",
            .publication = "10.1080/10629360701844142",
            .notes = "Requires a connected molecular bond graph because "
                     "every atom pair needs a finite bond-path distance.",
            .priority = 55};

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
