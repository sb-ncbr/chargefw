#include "methods/linear_algebra/blocked_factorization.h"
#include "methods/linear_algebra/dense_solve.h"

#include <Eigen/LU>

#include <algorithm>
#include <utility>

namespace chargefw::methods::linear_algebra {
namespace {

using detail::block_size;
using detail::blocked_size;
using detail::for_each_chunk;
using Eigen::Index;

using Transpositions = Eigen::Transpositions<Eigen::Dynamic, Eigen::Dynamic, int>;

// Partial-pivoting LU of a panel with at least as many rows as columns, by halving its columns:
// factorize the left half, bring the right half up to date, factorize what remains of it. Step i
// swapped row i with row transpositions[i] of the panel. A zero pivot column is left as it is, as
// in Eigen.
// NOLINTNEXTLINE(misc-no-recursion)
auto factorize_panel(Eigen::Ref<Eigen::MatrixXd> panel, int* transpositions) -> void {
    const auto rows = panel.rows();
    const auto columns = panel.cols();
    if (columns == 1) {
        Index pivot = 0;
        const auto largest = panel.col(0).cwiseAbs().maxCoeff(&pivot);
        *transpositions = static_cast<int>(pivot);
        if (largest != 0.0) {
            std::swap(panel(0, 0), panel(pivot, 0));
            panel.col(0).tail(rows - 1) /= panel(0, 0);
        }
        return;
    }

    const auto left = columns / 2;
    const auto right = columns - left;
    factorize_panel(panel.leftCols(left), transpositions);

    auto right_columns = panel.rightCols(right);
    for (Index row = 0; row < left; ++row) {
        right_columns.row(row).swap(right_columns.row(transpositions[row]));
    }
    auto top_right = panel.topRightCorner(left, right);
    panel.topLeftCorner(left, left).triangularView<Eigen::UnitLower>().solveInPlace(top_right);
    panel.bottomRightCorner(rows - left, right).noalias() -=
        panel.bottomLeftCorner(rows - left, left) * top_right;

    factorize_panel(panel.bottomRightCorner(rows - left, right), transpositions + left);

    auto left_columns = panel.leftCols(left);
    for (Index row = left; row < columns; ++row) {
        transpositions[row] += static_cast<int>(left);
        left_columns.row(row).swap(left_columns.row(transpositions[row]));
    }
}

// Right-looking blocked LU: factorize a panel of columns, apply its row swaps to the other columns,
// and bring the columns to its right up to date. Returns the row swapped with each row in turn.
[[nodiscard]] auto blocked_lu_in_place(Eigen::MatrixXd& matrix) -> Transpositions {
    const auto size = matrix.rows();
    auto transpositions = Transpositions{size};
    auto& swapped_with = transpositions.indices();
    for (Index start = 0; start < size; start += block_size) {
        const auto width = std::min(block_size, size - start);
        factorize_panel(matrix.block(start, start, size - start, width), &swapped_with(start));
        for (Index row = start; row < start + width; ++row) {
            swapped_with(row) += static_cast<int>(start);
        }

        const auto swap_rows = [&](const Index first, const Index columns) {
            auto chunk = matrix.middleCols(first, columns);
            for (Index row = start; row < start + width; ++row) {
                chunk.row(row).swap(chunk.row(swapped_with(row)));
            }
        };
        for_each_chunk(start, swap_rows);

        const auto after = start + width;
        const auto lower =
            matrix.block(start, start, width, width).triangularView<Eigen::UnitLower>();
        for_each_chunk(size - after, [&](const Index first, const Index columns) {
            swap_rows(after + first, columns);
            auto top = matrix.block(start, after + first, width, columns);
            lower.solveInPlace(top);
            matrix.block(after, after + first, size - after, columns).noalias() -=
                matrix.block(after, start, size - after, width) * top;
        });
    }
    return transpositions;
}

template <typename Rhs>
[[nodiscard]] auto solve_factored(const Eigen::MatrixXd& factors, Rhs x) -> Rhs {
    factors.triangularView<Eigen::UnitLower>().solveInPlace(x);
    factors.triangularView<Eigen::Upper>().solveInPlace(x);
    return x;
}

} // namespace

InPlaceLu::InPlaceLu(Eigen::MatrixXd& matrix) : factors_{&matrix} {
    if (matrix.rows() >= blocked_size) {
        row_permutation_ = blocked_lu_in_place(matrix);
    } else {
        row_permutation_ = Eigen::PartialPivLU<Eigen::Ref<Eigen::MatrixXd>>{matrix}.permutationP();
    }
}

auto InPlaceLu::solve(const Eigen::VectorXd& rhs) const -> Eigen::VectorXd {
    return solve_factored<Eigen::VectorXd>(*factors_, row_permutation_ * rhs);
}

auto InPlaceLu::solve(const Eigen::MatrixXd& rhs) const -> Eigen::MatrixXd {
    return solve_factored<Eigen::MatrixXd>(*factors_, row_permutation_ * rhs);
}

auto solve_in_place(Eigen::MatrixXd& matrix, const Eigen::VectorXd& rhs) -> Eigen::VectorXd {
    return InPlaceLu{matrix}.solve(rhs);
}

} // namespace chargefw::methods::linear_algebra
