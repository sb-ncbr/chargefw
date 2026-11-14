#include "methods/linear_algebra/dense_product.h"

namespace chargefw::methods::linear_algebra {

auto multiply(const Eigen::MatrixXd& left, const Eigen::MatrixXd& right) -> Eigen::MatrixXd {
    return left * right;
}

auto multiply_transposed(const Eigen::MatrixXd& left, const Eigen::MatrixXd& right)
    -> Eigen::MatrixXd {
    return left.transpose() * right;
}

} // namespace chargefw::methods::linear_algebra
