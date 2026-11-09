#pragma once

#include <cstddef>
#include <limits>
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

// Groups a method's parameter values into columns for LinearizedCharges::jacobian(). Each atom
// parameter value, bond parameter value, and common parameter belongs to one column, or to none
// with `none`. In a fit, a column is one value of a parameter set, shared by all atoms or bonds of
// a type. The parameters are those of Method::requirements(), in its order.
struct ParameterColumns {
    static constexpr std::size_t none = std::numeric_limits<std::size_t>::max();
    std::size_t count = 0;
    std::vector<std::vector<std::size_t>> atom; // [atom parameter][atom]
    std::vector<std::vector<std::size_t>> bond; // [bond parameter][bond]
    std::vector<std::size_t> common;            // [common parameter]
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

    // The derivatives of every charge with respect to grouped parameter values: column c is the
    // change of the charges when every value in column c changes by the same amount. Returns
    // atoms x columns.count values, column after column: entry [column * atom_count + atom]. Costs
    // about one solve per column, so it suits few columns, for example the parameter-set values of
    // one molecule. Throws std::invalid_argument unless `columns` has one entry per parameter value
    // and every column is below columns.count or none, and std::runtime_error when a value is not
    // finite. May be called concurrently with itself and with gradient().
    [[nodiscard]] virtual auto jacobian(const ParameterColumns& columns) const
        -> std::vector<double> = 0;

  protected:
    LinearizedCharges() = default;
};

} // namespace chargefw::methods
