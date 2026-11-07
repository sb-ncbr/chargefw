#pragma once

#include <span>
#include <vector>

namespace chargefw::methods {

// A gradient with respect to a method's parameter values, before they are aggregated into
// parameter-set entries: one value per atom for each atom parameter, one per bond for each bond
// parameter, and one per common parameter. The parameters are those of Method::requirements(), in
// its order.
struct ParameterGradient {
    std::vector<std::vector<double>> atom; // [atom parameter][atom]
    std::vector<std::vector<double>> bond; // [bond parameter][bond]
    std::vector<double> common;            // [common parameter]
};

// A method's charges at one input, with their exact first-order dependence on its parameter
// values. Holds the factorized system, so a gradient needs one more solve with the stored factor
// instead of a new calculation.
//
// The object owns everything it needs and stays valid after the input, its parameter set, and their
// storage are destroyed. gradient() may be called concurrently on one object. Geometry, topology,
// classification, target charge, and fixed sources are held constant.
class LinearizedCharges {
  public:
    LinearizedCharges(const LinearizedCharges&) = delete;
    auto operator=(const LinearizedCharges&) -> LinearizedCharges& = delete;
    LinearizedCharges(LinearizedCharges&&) = delete;
    auto operator=(LinearizedCharges&&) -> LinearizedCharges& = delete;
    virtual ~LinearizedCharges() = default;

    // The charges at the point of linearization, equal to Method::calculate() on the same input.
    // The span stays valid until the object is destroyed.
    [[nodiscard]] virtual auto charges() const noexcept -> std::span<const double> = 0;

    // The gradient of sum(cotangent[i] * charges[i]) with respect to every parameter value. With a
    // cotangent that is 1 at one atom and 0 elsewhere, it is that atom charge's gradient. Throws
    // std::invalid_argument unless `cotangent` has one finite value per atom, and
    // std::runtime_error when the gradient is not finite.
    [[nodiscard]] virtual auto gradient(std::span<const double> cotangent) const
        -> ParameterGradient = 0;

  protected:
    LinearizedCharges() = default;
};

} // namespace chargefw::methods
