#pragma once

#include <chargefw/features/prepared_molecule.h>
#include <chargefw/methods/method.h>

#include <array>

namespace chargefw::methods::builtin {

class Charge2Method final : public Method {
  public:
    [[nodiscard]] auto metadata() const noexcept -> const MethodMetadata& override {
        static constexpr MethodMetadata metadata{.id = "charge2",
                                                 .name = "Charge2",
                                                 .full_name = "Charge2",
                                                 .publication = "10.1002/jcc.540030316",
                                                 .priority = 30,
                                                 .recommended = false};

        return metadata;
    }

    [[nodiscard]] auto requirements() const -> MethodRequirements override {
        auto requirements = MethodRequirements{};
        requirements.bonds = true;
        requirements.common_parameters = {"a1", "a2", "a3", "b", "c", "alpha"};
        requirements.atom_parameters = {"chi", "P0", "q0"};
        requirements.resources.time = ComplexityTerm::atoms_plus_bonds;
        requirements.resources.memory = ComplexityTerm::atoms;
        return requirements;
    }

    [[nodiscard]] auto option_schema() const noexcept
        -> std::span<const MethodOptionSpec> override {
        static const std::array option_schema{
            MethodOptionSpec{.id = "iters",
                             .description = "Maximum number of Charge2 iterations",
                             .type = MethodOptionType::integer,
                             .default_value = 10,
                             .choices = {},
                             .minimum = 1,
                             .maximum = 25}};
        return {option_schema.data(), option_schema.size()};
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
