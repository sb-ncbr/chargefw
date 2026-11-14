#include "methods/linear_algebra/dense_solve.h"

#include <Eigen/Cholesky>

namespace chargefw::methods::linear_algebra {

auto cholesky_in_place(Eigen::MatrixXd& matrix) -> bool {
    const Eigen::LLT<Eigen::Ref<Eigen::MatrixXd>> llt{matrix};
    return llt.info() == Eigen::Success;
}

auto cholesky_solve_in_place(const Eigen::MatrixXd& factor, Eigen::VectorXd& rhs) -> void {
    const auto lower = factor.triangularView<Eigen::Lower>();
    lower.solveInPlace(rhs);
    lower.transpose().solveInPlace(rhs);
}

auto cholesky_solve_in_place(const Eigen::MatrixXd& factor, Eigen::MatrixXd& rhs) -> void {
    const auto lower = factor.triangularView<Eigen::Lower>();
    lower.solveInPlace(rhs);
    lower.transpose().solveInPlace(rhs);
}

} // namespace chargefw::methods::linear_algebra
