#include "methods/linear_algebra/dense_solve.h"

#include <Eigen/Cholesky>

#include <algorithm>

#include <oneapi/tbb/parallel_for.h>
#include <oneapi/tbb/task_arena.h>

namespace chargefw::methods::linear_algebra {
namespace {

using Eigen::Index;

// Matrices of at least this size are factorized block by block, with the threads of the current
// arena sharing the work of a block. Smaller ones, such as the fragments of cutoff and cover, go to
// Eigen's factorization, where threads would gain nothing.
constexpr Index blocked_size = 1500;
// Columns factorized per step, and rows or columns of the remaining matrix that one task updates.
// Tasks are cut by these sizes only, never by the number of threads, so the factor is the same for
// any number of threads.
constexpr Index block_size = 256;
constexpr Index chunk_size = 128;

// function(first, size) for consecutive chunks of [0, size), in parallel. Isolated, so that a
// thread waiting for its chunks does not start another molecule of the enclosing loop and allocate
// that molecule's matrices on top of this one.
template <typename Function>
auto for_each_chunk(const Index size, const Function& function) -> void {
    const auto chunk_count = (size + chunk_size - 1) / chunk_size;
    oneapi::tbb::this_task_arena::isolate([&] {
        oneapi::tbb::parallel_for(Index{0}, chunk_count, [&](const Index chunk) {
            const auto first = chunk * chunk_size;
            function(first, std::min(chunk_size, size - first));
        });
    });
}

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
