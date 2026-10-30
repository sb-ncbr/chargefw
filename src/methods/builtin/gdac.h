#pragma once

#include <chargefw/methods/calculation_input.h>
#include <chargefw/methods/method.h>
#include <chargefw/methods/method_metadata.h>
#include <chargefw/methods/method_requirements.h>

#include <array>

namespace chargefw::methods::builtin {

class GDACMethod final : public Method {
  public:
    [[nodiscard]] auto metadata() const noexcept -> const MethodMetadata& override {
        static constexpr MethodMetadata metadata{
            .id = "gdac",
            .name = "GDAC",
            .full_name = "Geometry-Dependent Net Atomic Charges",
            .publication = "10.1021/jp0023213",
            .notes = "This implementation omits charged-species "
                     "initialization; every connected component must be "
                     "formally neutral.",
            .priority = 100};

        return metadata;
    }

    [[nodiscard]] auto requirements() const -> MethodRequirements override {
        auto requirements = MethodRequirements{};
        requirements.coordinates = true;
        requirements.bonds = true;
        requirements.atom_parameters = {"A", "B"};
        requirements.resources.time = ComplexityTerm::atoms_plus_bonds;
        requirements.resources.memory = ComplexityTerm::atoms;
        return requirements;
    }

    [[nodiscard]] auto option_schema() const noexcept
        -> std::span<const MethodOptionSpec> override {
        static const std::array option_schema{
            MethodOptionSpec{.id = "iters",
                             .description = "Number of GDAC iterations",
                             .type = MethodOptionType::integer,
                             .default_value = 7,
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
