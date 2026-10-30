#pragma once

#include <chargefw/methods/method.h>

namespace chargefw::methods::builtin {

class SQEq0Method final : public Method {
  public:
    [[nodiscard]] auto metadata() const noexcept -> const MethodMetadata& override {
        static constexpr MethodMetadata metadata{
            .id = "sqeq0",
            .name = "SQE+q0",
            .full_name = "Split-charge Equilibration with Initial Formal Charges",
            .publication = "10.1021/ct200512e",
            .notes = "In full execution, atomic formal charges define the reference state and each "
                     "connected component preserves its formal-charge total.",
            .priority = 80};

        return metadata;
    }

    [[nodiscard]] auto requirements() const -> MethodRequirements override {
        auto requirements = MethodRequirements{};
        requirements.coordinates = true;
        requirements.bonds = true;
        requirements.atom_parameters = {"electronegativity", "hardness", "width"};
        requirements.bond_parameters = {"kappa"};
        requirements.resources.time = ComplexityTerm::bonds_cubed;
        requirements.resources.memory = ComplexityTerm::atoms_squared_plus_bonds_squared;
        requirements.resources.supports_cutoff = true;
        requirements.resources.supports_cover = true;
        requirements.resources.reduced_charge_policy =
            ReducedChargePolicy::formal_charge_components;
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
