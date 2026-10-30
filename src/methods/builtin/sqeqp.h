#pragma once

#include <chargefw/methods/method.h>

namespace chargefw::methods::builtin {

class SQEqpMethod final : public Method {
  public:
    [[nodiscard]] auto metadata() const noexcept -> const MethodMetadata& override {
        static constexpr MethodMetadata metadata{
            .id = "sqeqp",
            .name = "SQE+qp",
            .full_name = "Split-charge Equilibration with Parameterized Initial Charges",
            .publication = "10.1186/s13321-021-00528-w",
            .notes = "Parameterized reference charges are normalized globally to the calculation "
                     "target charge before equilibration.",
            .priority = 210};

        return metadata;
    }

    [[nodiscard]] auto requirements() const -> MethodRequirements override {
        auto requirements = MethodRequirements{};
        requirements.coordinates = true;
        requirements.bonds = true;
        requirements.atom_parameters = {"electronegativity", "hardness", "width", "q0"};
        requirements.bond_parameters = {"kappa"};
        requirements.resources.time = ComplexityTerm::bonds_cubed;
        requirements.resources.memory = ComplexityTerm::atoms_squared_plus_bonds_squared;
        requirements.resources.supports_cutoff = true;
        requirements.resources.supports_cover = true;
        requirements.resources.reduced_charge_policy =
            ReducedChargePolicy::parameterized_charge_components;
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
