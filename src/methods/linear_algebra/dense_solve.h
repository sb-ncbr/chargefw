#pragma once

#include <Eigen/Core>
#include <Eigen/LU>

namespace chargefw::methods::linear_algebra {

// The dense factorizations are defined in dense_lu.cpp and dense_cholesky.cpp, so that Eigen's
// blocked LU and Cholesky code, with the matrix-product kernels they use, is compiled in those two
// files instead of in every method that solves a system.

// Partial-pivoting LU factors of a square matrix. The matrix is overwritten by its factors instead
// of copied, which halves peak memory and cache footprint for large dense systems; it must outlive
// the factorization.
class InPlaceLu {
  public:
    explicit InPlaceLu(Eigen::MatrixXd& matrix);

    [[nodiscard]] auto solve(const Eigen::VectorXd& rhs) const -> Eigen::VectorXd;
    [[nodiscard]] auto solve(const Eigen::MatrixXd& rhs) const -> Eigen::MatrixXd;

  private:
    Eigen::PartialPivLU<Eigen::Ref<Eigen::MatrixXd>> lu_;
};

// Solves matrix * x = rhs with partial-pivoting LU, overwriting the matrix by its factors.
[[nodiscard]] auto solve_in_place(Eigen::MatrixXd& matrix, const Eigen::VectorXd& rhs)
    -> Eigen::VectorXd;

// Replaces the lower triangle of a symmetric matrix by its Cholesky factor L, reading only that
// triangle. Returns false when the matrix is not positive definite to working precision.
[[nodiscard]] auto cholesky_in_place(Eigen::MatrixXd& matrix) -> bool;

// Solves L L^T x = rhs in place, for the factor that cholesky_in_place() left in the lower
// triangle.
auto cholesky_solve_in_place(const Eigen::MatrixXd& factor, Eigen::VectorXd& rhs) -> void;
auto cholesky_solve_in_place(const Eigen::MatrixXd& factor, Eigen::MatrixXd& rhs) -> void;

} // namespace chargefw::methods::linear_algebra
