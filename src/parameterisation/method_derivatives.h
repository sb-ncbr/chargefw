#pragma once

#include <chargefw/methods/calculation_input.h>
#include <chargefw/methods/method.h>

#include <cstddef>
#include <limits>
#include <memory>
#include <span>
#include <vector>

// The exact derivatives of methods' charges with respect to their parameter values, private to the
// library. They are attached to the built-in methods from outside: the methods know nothing of
// them, and share with them only the systems they solve (sqe_system.h, eem_system.h).
namespace chargefw::parameterisation::detail {

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
// values. It keeps what the calculation produced, so that derivatives cost less than calculating
// the charges again.
//
// The object owns everything it needs and stays valid after the input, its parameter set, and their
// storage are destroyed. jacobian() may be called concurrently on one object. Geometry, topology,
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

    // The derivatives of every charge with respect to grouped parameter values: column c is the
    // change of the charges per unit change of every value in column c together. Returns
    // atoms x columns.count values, column after column: entry [column * atom_count + atom]. Its
    // cost grows with the number of columns, so it suits few of them, such as the parameter-set
    // values of one molecule. Throws std::invalid_argument unless `columns` has one entry per
    // parameter value and every column is below columns.count or none, and std::runtime_error when
    // a value is not finite.
    [[nodiscard]] virtual auto jacobian(const ParameterColumns& columns) const
        -> std::vector<double> = 0;

  protected:
    LinearizedCharges() = default;
};

// The exact derivatives of one method's charges.
class MethodDerivatives {
  public:
    MethodDerivatives(const MethodDerivatives&) = delete;
    auto operator=(const MethodDerivatives&) -> MethodDerivatives& = delete;
    MethodDerivatives(MethodDerivatives&&) = delete;
    auto operator=(MethodDerivatives&&) -> MethodDerivatives& = delete;
    virtual ~MethodDerivatives() = default;

    // The charges at `input`, equal to Method::calculate() of the method, with their derivatives.
    // Throws for invalid input and numerical failures as calculate() does.
    [[nodiscard]] virtual auto linearize(const methods::CalculationInput& input) const
        -> std::unique_ptr<LinearizedCharges> = 0;

  protected:
    MethodDerivatives() = default;
};

// The derivatives of `method`, or nullptr when it has none. Matched by the method's concrete
// built-in type, not by its ID.
[[nodiscard]] auto derivatives_of(const methods::Method& method) -> const MethodDerivatives*;

// The derivatives of each built-in method that has them (eem_derivatives.cpp,
// sqe_derivatives.cpp).
[[nodiscard]] auto eem_derivatives() -> const MethodDerivatives&;
[[nodiscard]] auto sqe_derivatives() -> const MethodDerivatives&;
[[nodiscard]] auto sqeq0_derivatives() -> const MethodDerivatives&;
[[nodiscard]] auto sqeqp_derivatives() -> const MethodDerivatives&;

} // namespace chargefw::parameterisation::detail
