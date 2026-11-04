#pragma once

#include <chargefw/methods/method.h>

#include <array>

namespace chargefw::methods::builtin {

class MPEOEMethod final : public Method {
  public:
    [[nodiscard]] auto metadata() const noexcept -> const MethodMetadata& override {
        static constexpr MethodMetadata metadata{
            .id = "mpeoe",
            .name = "MPEOE",
            .full_name = "Modified Partial Equalization of Orbital Electronegativity",
            .publication = "10.1021/j100374a066",
            .notes = "With the default initial_charges=zero, the formal charges of every bonded "
                     "part must add up to zero. initial_charges=formal, added by ChargeFW, starts "
                     "from the formal charges and allows charged molecules.",
            .priority = 110};

        return metadata;
    }

    [[nodiscard]] auto requirements() const -> MethodRequirements override {
        auto requirements = MethodRequirements{};
        requirements.bonds = true;
        requirements.common_parameters = {"Hplus"};
        requirements.atom_parameters = {"A", "B"};
        requirements.bond_parameters = {"f"};
        requirements.resources.time = ComplexityTerm::atoms_plus_bonds;
        requirements.resources.memory = ComplexityTerm::atoms;
        return requirements;
    }

    [[nodiscard]] auto option_schema() const noexcept
        -> std::span<const MethodOptionSpec> override {

        static const std::array option_schema{
            MethodOptionSpec{.id = "iters",
                             .description = "Number of MPEOE iterations",
                             .type = MethodOptionType::integer,
                             .default_value = 7,
                             .choices = {},
                             .minimum = 1,
                             .maximum = 25},
            MethodOptionSpec{.id = "initial_charges",
                             .description = "Atomic charges before MPEOE charge transfer",
                             .type = MethodOptionType::string,
                             .default_value = std::string{"zero"},
                             .choices = {std::string{"zero"}, std::string{"formal"}}}};
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
