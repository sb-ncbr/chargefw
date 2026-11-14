#include "methods/linear_algebra/dense_solve.h"

namespace chargefw::methods::linear_algebra {

InPlaceLu::InPlaceLu(Eigen::MatrixXd& matrix) : lu_{matrix} {}

auto InPlaceLu::solve(const Eigen::VectorXd& rhs) const -> Eigen::VectorXd {
    return lu_.solve(rhs);
}

auto InPlaceLu::solve(const Eigen::MatrixXd& rhs) const -> Eigen::MatrixXd {
    return lu_.solve(rhs);
}

auto solve_in_place(Eigen::MatrixXd& matrix, const Eigen::VectorXd& rhs) -> Eigen::VectorXd {
    return InPlaceLu{matrix}.solve(rhs);
}

} // namespace chargefw::methods::linear_algebra
