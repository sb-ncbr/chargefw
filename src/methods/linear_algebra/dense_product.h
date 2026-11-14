#pragma once

#include <Eigen/Core>

namespace chargefw::methods::linear_algebra {

// Products of two dense matrices, defined in dense_product.cpp so that Eigen's matrix-product
// kernels are compiled there instead of in every file that multiplies matrices. The arguments are
// evaluated into matrices first, so pass expressions only where that temporary is acceptable.

// left * right.
[[nodiscard]] auto multiply(const Eigen::MatrixXd& left, const Eigen::MatrixXd& right)
    -> Eigen::MatrixXd;
// left^T * right.
[[nodiscard]] auto multiply_transposed(const Eigen::MatrixXd& left, const Eigen::MatrixXd& right)
    -> Eigen::MatrixXd;

} // namespace chargefw::methods::linear_algebra
