#include "methods/linear_algebra/dense_solve.h"

#include <Eigen/Cholesky>
#include <Eigen/Core>

#include <snitch/snitch.hpp>

#include <oneapi/tbb/task_arena.h>

namespace linear_algebra = chargefw::methods::linear_algebra;

namespace {

// Large enough for the blocked factorization, with a last block and a last chunk that are not full.
constexpr Eigen::Index blocked_size = 1600;

// Symmetric and diagonally dominant, so positive definite.
auto make_positive_definite(const Eigen::Index size) -> Eigen::MatrixXd {
    Eigen::MatrixXd lower = Eigen::MatrixXd::Zero(size, size);
    for (Eigen::Index column = 0; column < size; ++column) {
        lower(column, column) = static_cast<double>(size);
        for (Eigen::Index row = column + 1; row < size; ++row) {
            lower(row, column) = 1.0 / static_cast<double>(1 + row - column + (row * column) % 7);
        }
    }
    return lower.selfadjointView<Eigen::Lower>();
}

auto factorize(Eigen::MatrixXd matrix, const int threads) -> Eigen::MatrixXd {
    oneapi::tbb::task_arena arena{threads};
    const auto positive_definite =
        arena.execute([&] { return linear_algebra::cholesky_in_place(matrix); });
    REQUIRE(positive_definite);
    return matrix;
}

} // namespace

TEST_CASE("the Cholesky factor of a large matrix solves its system for any number of threads",
          "[methods][linear-algebra]") {
    const auto matrix = make_positive_definite(blocked_size);
    const Eigen::VectorXd rhs = Eigen::VectorXd::LinSpaced(blocked_size, -1.0, 1.0);

    const auto factor = factorize(matrix, 1);
    Eigen::VectorXd solution = rhs;
    linear_algebra::cholesky_solve_in_place(factor, solution);
    const Eigen::VectorXd expected = matrix.llt().solve(rhs);
    CHECK((solution - expected).cwiseAbs().maxCoeff() < 1e-12);

    // The upper triangle is left as it was.
    const Eigen::MatrixXd upper = factor.triangularView<Eigen::StrictlyUpper>();
    const Eigen::MatrixXd original_upper = matrix.triangularView<Eigen::StrictlyUpper>();
    CHECK(upper == original_upper);

    CHECK(factorize(matrix, 4) == factor);
}

TEST_CASE("a large matrix that is not positive definite has no Cholesky factor",
          "[methods][linear-algebra]") {
    auto matrix = make_positive_definite(blocked_size);
    // In the last block, after full blocks were factorized.
    matrix(blocked_size - 2, blocked_size - 2) = -1.0;

    CHECK_FALSE(linear_algebra::cholesky_in_place(matrix));
}
