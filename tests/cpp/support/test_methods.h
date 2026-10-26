#pragma once

#include <chargefw/charges/atomic_charges.h>
#include <chargefw/methods/calculation_input.h>
#include <chargefw/methods/method.h>
#include <chargefw/methods/method_metadata.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/methods/method_requirements.h>

#include <optional>
#include <span>
#include <string_view>
#include <utility>
#include <vector>

namespace chargefw::test {

/// Test method with fixed metadata and requirements, no options, and zero charges. Derived test
/// methods override calculate() when a scenario needs other charges or side effects. The ID must
/// outlive the method, which string literals do.
class StubMethod : public methods::Method {
  public:
    explicit StubMethod(const std::string_view id, methods::MethodRequirements requirements = {},
                        const int priority = 0)
        : metadata_{.id = id,
                    .name = id,
                    .full_name = id,
                    .publication = std::nullopt,
                    .priority = priority},
          requirements_{std::move(requirements)} {}

    [[nodiscard]] auto metadata() const noexcept -> const methods::MethodMetadata& override {
        return metadata_;
    }

    [[nodiscard]] auto requirements() const -> methods::MethodRequirements override {
        return requirements_;
    }

    [[nodiscard]] auto option_schema() const noexcept
        -> std::span<const methods::MethodOptionSpec> override {
        return {};
    }

    [[nodiscard]] auto calculate(const methods::CalculationInput& input) const
        -> charges::AtomicCharges override {
        return charges::AtomicCharges{std::vector<double>(input.molecule().atom_count(), 0.0)};
    }

    [[nodiscard]] auto estimate(const methods::TargetSize& target,
                                const methods::MethodOptions& /*options*/) const
        -> methods::ResourceEstimate override {
        return {.memory_bytes = 0.0, .operations = target.atoms};
    }

  private:
    methods::MethodMetadata metadata_;
    methods::MethodRequirements requirements_;
};

} // namespace chargefw::test
