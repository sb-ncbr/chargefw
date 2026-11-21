#pragma once

#include <cstddef>

#include <oneapi/tbb/parallel_for.h>
#include <oneapi/tbb/task_arena.h>

namespace chargefw::methods::detail {

// From this many rows, the rows of a matrix are filled in parallel. Below, as for the fragments of
// cutoff and cover, threads would gain nothing.
inline constexpr std::size_t parallel_row_count = 1500;

// function(row) for each row in [0, count), on the threads of the current arena when there are
// many rows. The function must write only what belongs to its row. Isolated, so that a waiting
// thread does not start another molecule of the enclosing loop.
template <typename Function>
auto for_each_row(const std::size_t count, const Function& function) -> void {
    if (count < parallel_row_count) {
        for (std::size_t row = 0; row < count; ++row) {
            function(row);
        }
        return;
    }
    oneapi::tbb::this_task_arena::isolate(
        [&] { oneapi::tbb::parallel_for(std::size_t{0}, count, function); });
}

} // namespace chargefw::methods::detail
