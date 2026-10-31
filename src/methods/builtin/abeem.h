#pragma once

#include <chargefw/methods/method.h>

namespace chargefw::methods::builtin {

class ABEEMMethod final : public Method {
  public:
    [[nodiscard]] auto metadata() const noexcept -> const MethodMetadata& override {
        static constexpr MethodMetadata metadata{
            .id = "abeem",
            .name = "ABEEM",
            .full_name = "Atom-Bond Electronegativity Equalization Method",
            .publication = "10.1021/jp9711048",
            .priority = 190};

        return metadata;
    }

    [[nodiscard]] auto requirements() const -> MethodRequirements override {
        auto requirements = MethodRequirements{};
        requirements.coordinates = true;
        requirements.bonds = true;
        requirements.common_parameters = {"k"};
        requirements.atom_parameters = {"a", "b", "c"};
        requirements.bond_parameters = {"A", "B", "C", "D"};
        requirements.resources.time = ComplexityTerm::atoms_plus_bonds_cubed;
        requirements.resources.memory = ComplexityTerm::atoms_plus_bonds_squared;
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

  protected:
    auto add_method_specific_prerequisite_issues(const MethodPrerequisiteInput& input,
                                                 PrerequisiteResult& result) const -> void override;
};

} // namespace chargefw::methods::builtin
