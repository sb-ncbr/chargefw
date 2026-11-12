#pragma once

#include "calculation/calculation_cancelled.h"
#include "calculation/observer_notifications.h"

#include <chargefw/calculation/observer.h>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <mutex>
#include <optional>
#include <string_view>
#include <utility>

#include <oneapi/tbb/global_control.h>
#include <oneapi/tbb/parallel_for.h>
#include <oneapi/tbb/task_arena.h>

namespace chargefw::calculation::detail {

// A thread limit caps oneTBB's allowed parallelism, which bounds the workers any arena can use.
// Zero uses that parallelism.
[[nodiscard]] inline auto thread_count(const std::size_t max_threads) -> std::size_t {
    const auto allowed = oneapi::tbb::global_control::active_value(
        oneapi::tbb::global_control::max_allowed_parallelism);
    return max_threads == 0 ? allowed : std::min(max_threads, allowed);
}

// The callable must write only to the result associated with its index; execution order is
// intentionally unspecified.
template <typename Function>
auto parallel_for_indexed(const std::size_t count, const std::size_t max_threads,
                          Function&& function) -> void {
    if (count == 0) {
        return;
    }

    const auto threads = thread_count(max_threads);
    if (threads == 1) {
        for (std::size_t index = 0; index < count; ++index) {
            function(index);
        }
        return;
    }

    oneapi::tbb::task_arena arena{static_cast<int>(threads)};
    arena.execute([&] {
        oneapi::tbb::parallel_for(std::size_t{0}, count,
                                  [&function](const std::size_t index) { function(index); });
    });
}

// Context carried through a parallel loop for progress emission. observer, mode, and method_id are
// shared; target_index/molecule_index/conformer_index identify the enclosing target.
struct ProgressContext {
    // NOLINTNEXTLINE(cppcoreguidelines-avoid-const-or-ref-data-members)
    const CalculationObserver& observer;
    ExecutionMode mode{};
    std::string_view method_id;
    std::size_t target_index{};
    std::size_t target_count{};
    std::size_t molecule_index{};
    std::optional<std::size_t> conformer_index{};
    std::chrono::steady_clock::time_point computation_start{};
};

// Throws CalculationCancelled if the observer requests termination. Called by executors at
// cancellation check points (start of each target iteration).
inline auto check_cancellation(const CalculationObserver& observer) -> void {
    if (observer.cancelled()) {
        throw CalculationCancelled{};
    }
}

// Emits a single target-tier progress event. Called by executors at the start and end of each
// target iteration, where the per-target molecule_index and conformer_index are known.
inline auto emit_target_event(const CalculationPhase phase, const ProgressContext& ctx) -> void {
    report_progress(
        ctx.observer,
        CalculationProgress{
            .phase = phase,
            .mode = ctx.mode,
            .method_id = ctx.method_id,
            .target_index = ctx.target_index,
            .target_count = ctx.target_count,
            .molecule_index = ctx.molecule_index,
            .conformer_index = ctx.conformer_index,
            .elapsed_seconds = std::chrono::duration<double>{std::chrono::steady_clock::now() -
                                                             ctx.computation_start}
                                   .count(),
        });
}

// Minimum interval between fragment-tier progress events, in nanoseconds.
inline constexpr std::int64_t fragment_throttle_ns = 200'000'000; // 200 ms

// Fragment-tier loop with throttled progress emission and cooperative cancellation. Checks
// cancellation at the start of each iteration. The first completed fragment is reported
// immediately; subsequent snapshots are throttled and a terminal snapshot is emitted only when
// the throttled snapshots have not already reported full completion. Snapshots for an individual
// target are serialized, so their completed_fragment_count values strictly increase. Callbacks for
// different targets can run concurrently, so observers must not use global callback arrival order
// to infer completion order.
//
// Target-tier events are NOT emitted here — executors emit them via check_cancellation() and
// emit_target_event() so per-target molecule/conformer identity is correctly populated.
template <typename Function>
auto progress_for_indexed(const std::size_t count, const std::size_t max_threads,
                          const ProgressContext& ctx, Function&& function) -> void {
    if (count == 0) {
        parallel_for_indexed(count, max_threads, std::forward<Function>(function));
        return;
    }

    const auto& observer = ctx.observer;

    const auto elapsed = [&ctx]() -> double {
        return std::chrono::duration<double>{std::chrono::steady_clock::now() -
                                             ctx.computation_start}
            .count();
    };

    constexpr auto no_fragment_event = std::numeric_limits<std::int64_t>::min();
    std::atomic<std::int64_t> last_emit_ns{no_fragment_event};
    std::atomic<std::size_t> completed{0};
    auto emission_mutex = std::mutex{};
    auto last_emitted_completed = std::size_t{0};

    const auto emit_progress = [&](const std::size_t done) {
        const auto lock = std::scoped_lock{emission_mutex};
        if (done <= last_emitted_completed) {
            return;
        }
        last_emitted_completed = done;
        report_progress(observer, CalculationProgress{
                                      .phase = CalculationPhase::fragment_progress,
                                      .mode = ctx.mode,
                                      .method_id = ctx.method_id,
                                      .target_index = ctx.target_index,
                                      .target_count = ctx.target_count,
                                      .completed_fragment_count = done,
                                      .fragment_count = count,
                                      .molecule_index = ctx.molecule_index,
                                      .conformer_index = ctx.conformer_index,
                                      .elapsed_seconds = elapsed(),
                                  });
    };

    auto wrapped = [&](const std::size_t index) {
        if (observer.cancelled()) {
            throw CalculationCancelled{};
        }
        function(index);
        const auto done = completed.fetch_add(1, std::memory_order_relaxed) + 1;

        const auto now_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(
                                std::chrono::steady_clock::now().time_since_epoch())
                                .count();
        auto last = last_emit_ns.load(std::memory_order_relaxed);
        if ((last == no_fragment_event || now_ns - last >= fragment_throttle_ns) &&
            last_emit_ns.compare_exchange_strong(last, now_ns, std::memory_order_relaxed)) {
            emit_progress(done);
        }
    };
    parallel_for_indexed(count, max_threads, wrapped);

    if (last_emitted_completed != count) {
        emit_progress(count);
    }
}

} // namespace chargefw::calculation::detail
