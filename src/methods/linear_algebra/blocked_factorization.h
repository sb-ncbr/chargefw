#pragma once

#include <Eigen/Core>

#include <algorithm>

#include <oneapi/tbb/parallel_for.h>
#include <oneapi/tbb/task_arena.h>

namespace chargefw::methods::linear_algebra::detail {

// Matrices of at least this size are factorized block by block, with the threads of the current
// arena sharing the work of a block. Smaller ones, such as the fragments of cutoff and cover, go to
// Eigen's factorization, where threads would gain nothing.
inline constexpr Eigen::Index blocked_size = 1500;
// Columns factorized per step, and rows or columns of the remaining matrix that one task updates.
// Tasks are cut by these sizes only, never by the number of threads, so the factors are the same
// for any number of threads.
inline constexpr Eigen::Index block_size = 256;
inline constexpr Eigen::Index chunk_size = 128;

// function(first, size) for consecutive chunks of [0, size), in parallel. Isolated, so that a
// thread waiting for its chunks does not start another molecule of the enclosing loop and allocate
// that molecule's matrices on top of this one.
template <typename Function>
auto for_each_chunk(const Eigen::Index size, const Function& function) -> void {
    const auto chunk_count = (size + chunk_size - 1) / chunk_size;
    oneapi::tbb::this_task_arena::isolate([&] {
        oneapi::tbb::parallel_for(Eigen::Index{0}, chunk_count, [&](const Eigen::Index chunk) {
            const auto first = chunk * chunk_size;
            function(first, std::min(chunk_size, size - first));
        });
    });
}

} // namespace chargefw::methods::linear_algebra::detail
