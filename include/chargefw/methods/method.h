#pragma once

#include <chargefw/charges/atomic_charges.h>
#include <chargefw/methods/calculation_input.h>
#include <chargefw/methods/method_metadata.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/methods/method_prerequisites.h>
#include <chargefw/methods/method_requirements.h>

#include <span>
#include <string_view>

namespace chargefw::methods {

class Method {
  public:
    Method(const Method&) = delete;
    auto operator=(const Method&) -> Method& = delete;

    Method(Method&&) = delete;
    auto operator=(Method&&) -> Method& = delete;

    virtual ~Method() = default;

    [[nodiscard]] virtual auto metadata() const noexcept -> const MethodMetadata& = 0;

    [[nodiscard]] auto id() const noexcept -> std::string_view {
        return metadata().id;
    }

    [[nodiscard]] virtual auto requirements() const -> MethodRequirements = 0;

    [[nodiscard]] auto requires_parameters() const -> bool;

    [[nodiscard]] virtual auto option_schema() const noexcept
        -> std::span<const MethodOptionSpec> = 0;

    [[nodiscard]] auto check_method_prerequisites(const MethodPrerequisiteInput& input) const
        -> PrerequisiteResult;

    [[nodiscard]] virtual auto calculate(const CalculationInput& input) const
        -> charges::AtomicCharges = 0;

    // Options are validated and complete. The leading-order terms match requirements().resources.
    [[nodiscard]] virtual auto estimate(const TargetSize& target,
                                        const MethodOptions& options) const -> ResourceEstimate = 0;

  protected:
    Method() = default;

    virtual auto add_method_specific_prerequisite_issues(const MethodPrerequisiteInput& input,
                                                         PrerequisiteResult& result) const -> void;
};

} // namespace chargefw::methods
