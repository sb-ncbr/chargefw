#pragma once

#include <cstddef>
#include <limits>
#include <span>
#include <vector>

namespace chargefw::methods {

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

} // namespace chargefw::methods
