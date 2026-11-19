#pragma once

#include <cstdint>
#include <stdexcept>
#include <string>
#include <string_view>

namespace chargefw::methods {

// Why a method could not give charges for its input.
enum class FailureKind : std::uint8_t {
    // The calculation could not produce a stable, finite solution, for example because its
    // equations are singular or not positive definite.
    no_solution,
    // An iterative calculation did not converge.
    not_converged,
    // A parameter value is outside the range the method accepts.
    invalid_parameter_value,
};

[[nodiscard]] auto to_string(FailureKind value) -> std::string_view;

// What a method throws when it cannot give charges for its input. The message gives only the
// reason and what to do; the calculation adds the molecule, conformer, and fragment.
class MethodFailure : public std::runtime_error {
  public:
    MethodFailure(FailureKind kind, const std::string& message);

    [[nodiscard]] auto kind() const noexcept -> FailureKind;

  private:
    FailureKind kind_;
};

} // namespace chargefw::methods
