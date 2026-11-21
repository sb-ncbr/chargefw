#include "methods/linear_algebra/blocked_factorization.h"
#include "methods/linear_algebra/dense_solve.h"

#include <Eigen/Cholesky>

#include <algorithm>

namespace chargefw::methods::linear_algebra {
namespace {

using Eigen::Index;

using detail::block_size;
using detail::blocked_size;
using detail::for_each_chunk;

// Right-looking blocked Cholesky on the lower triangle: factorize a diagonal block, solve for the
// columns below it, and subtract their product from the rest of the lower triangle.
[[nodiscard]] auto blocked_cholesky_in_place(Eigen::MatrixXd& matrix) -> bool {
    const auto size = matrix.rows();
    for (Index start = 0; start < size; start += block_size) {
        const auto width = std::min(block_size, size - start);
        Eigen::Ref<Eigen::MatrixXd> diagonal = matrix.block(start, start, width, width);
        if (Eigen::LLT<Eigen::Ref<Eigen::MatrixXd>>{diagonal}.info() != Eigen::Success) {
            return false;
        }

        const auto below = start + width;
        const auto remaining = size - below;
        const auto factor = matrix.block(start, start, width, width).triangularView<Eigen::Lower>();
        for_each_chunk(remaining, [&](const Index first, const Index rows) {
            auto panel = matrix.block(below + first, start, rows, width);
            factor.transpose().solveInPlace<Eigen::OnTheRight>(panel);
        });
        for_each_chunk(remaining, [&](const Index first, const Index columns) {
            const auto corner = below + first;
            const auto lower = corner + columns;
            const auto panel = matrix.block(corner, start, columns, width);
            matrix.block(corner, corner, columns, columns)
                .selfadjointView<Eigen::Lower>()
                .rankUpdate(panel, -1.0);
            matrix.block(lower, corner, size - lower, columns).noalias() -=
                matrix.block(lower, start, size - lower, width) * panel.transpose();
        });
    }
    return true;
}

} // namespace

auto cholesky_in_place(Eigen::MatrixXd& matrix) -> bool {
    if (matrix.rows() >= blocked_size) {
        return blocked_cholesky_in_place(matrix);
    }
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
