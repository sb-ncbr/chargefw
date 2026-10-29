#include "cli_support.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <csignal>
#include <cstddef>
#include <exception>
#include <iostream>
#include <memory>
#include <mutex>
#include <optional>
#include <print>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace {

std::atomic_flag interruption_requested;

extern "C" void request_interruption(int /*signal*/) {
    interruption_requested.test_and_set(std::memory_order_relaxed);
}

class InterruptibleCalculationObserver final : public chargefw::calculation::CalculationObserver {
  public:
    explicit InterruptibleCalculationObserver(
        const chargefw::calculation::CalculationObserver* observer)
        : observer_{observer} {}

    void on_progress(const chargefw::calculation::CalculationProgress& progress) const override {
        if (observer_ != nullptr) {
            observer_->on_progress(progress);
        }
    }

    [[nodiscard]] auto cancelled() const noexcept -> bool override {
        return interruption_requested.test(std::memory_order_relaxed) ||
               (observer_ != nullptr && observer_->cancelled());
    }

  private:
    const chargefw::calculation::CalculationObserver* observer_;
};

class ScopedInterruptHandler final {
  public:
    ScopedInterruptHandler() {
        interruption_requested.clear(std::memory_order_relaxed);
        previous_handler_ = std::signal(SIGINT, request_interruption);
    }

    ScopedInterruptHandler(const ScopedInterruptHandler&) = delete;
    auto operator=(const ScopedInterruptHandler&) -> ScopedInterruptHandler& = delete;
    ScopedInterruptHandler(ScopedInterruptHandler&&) = delete;
    auto operator=(ScopedInterruptHandler&&) -> ScopedInterruptHandler& = delete;

    ~ScopedInterruptHandler() {
        std::signal(SIGINT, previous_handler_);
    }

  private:
    using SignalHandler = void (*)(int);

    SignalHandler previous_handler_ = SIG_DFL;
};

class TerminalProgressObserver final : public chargefw::calculation::CalculationObserver {
  public:
    void on_progress(const chargefw::calculation::CalculationProgress& progress) const override {
        const std::scoped_lock lock{mutex_};
        switch (progress.phase) {
        case chargefw::calculation::CalculationPhase::computation_started:
            std::print(std::cerr, "Calculating {}...\n", progress.method_id);
            break;
        case chargefw::calculation::CalculationPhase::target_finished:
            ++completed_targets_;
            render("Targets", completed_targets_, progress.target_count);
            break;
        case chargefw::calculation::CalculationPhase::fragment_progress:
            render("Fragments", progress.completed_fragment_count, progress.fragment_count,
                   progress.target_index + 1, progress.target_count);
            break;
        case chargefw::calculation::CalculationPhase::computation_finished:
            erase_line();
            std::print(std::cerr, "\n");
            break;
        default:
            break;
        }
        std::cerr.flush();
    }

  private:
    static constexpr std::size_t bar_width = 30;

    static void erase_line() {
        std::print(std::cerr, "\r\x1b[2K");
    }

    void render(const std::string_view label, const std::size_t complete, const std::size_t total,
                const std::optional<std::size_t> target_index = std::nullopt,
                const std::optional<std::size_t> target_count = std::nullopt) const {
        const auto filled = total == 0 ? 0U : std::min(bar_width, complete * bar_width / total);
        auto bar = std::string(bar_width, ' ');
        std::fill_n(bar.begin(), filled, '=');
        erase_line();
        if (target_index.has_value() && target_count.has_value()) {
            std::print(std::cerr, "{} {}/{} [{}] {}/{}", label, *target_index, *target_count, bar,
                       complete, total);
            return;
        }
        std::print(std::cerr, "{} [{}] {}/{}", label, bar, complete, total);
    }

    mutable std::mutex mutex_;
    mutable std::size_t completed_targets_ = 0;
};

auto run(std::span<char*> arguments) -> int {
    CLI::App app{"ChargeFW molecular charge calculation and inspection."};
    app.require_subcommand(1, 1);
    chargefw::cli::InputArguments calculate_input;
    chargefw::cli::InputArguments inspect_input;
    chargefw::cli::InputArguments applicability_input;
    chargefw::cli::SelectionArguments calculate_selection;
    chargefw::cli::SelectionArguments applicability_selection;
    chargefw::cli::OutputArguments calculate_output;
    std::string output_directory;
    std::string method_info;
    std::string parameter_method;
    std::string parameter_set_info;
    bool parameter_keys = false;
    bool progress = false;
    auto* calculate = app.add_subcommand("calculate", "Calculate and write partial charges");
    auto* inspect = app.add_subcommand("inspect", "Inspect imported molecular records");
    auto* applicability = app.add_subcommand("applicability", "Report applicable charge methods");
    auto* methods = app.add_subcommand("methods", "List registered charge methods");
    auto* parameters = app.add_subcommand("parameters", "List bundled parameter sets");
    chargefw::cli::add_input_options(*calculate, calculate_input);
    calculate->add_option("output", output_directory, "Output directory")->required();
    chargefw::cli::add_selection_options(*calculate, calculate_selection);
    calculate->add_option("--threads", calculate_selection.max_threads,
                          "Maximum calculation threads; 0 uses the oneTBB default");
    calculate->add_flag("--progress", progress, "Show calculation progress on standard error");
    calculate->add_flag("--output-mol2", calculate_output.mol2, "Write generated MOL2 output");
    calculate->add_flag("--output-mmcif", calculate_output.mmcif, "Write generated mmCIF output");
    calculate->add_flag("--output-pqr", calculate_output.pqr, "Write generated PQR output");
    chargefw::cli::add_input_options(*inspect, inspect_input);
    chargefw::cli::add_input_options(*applicability, applicability_input);
    chargefw::cli::add_selection_options(*applicability, applicability_selection);
    applicability->add_option("--threads", applicability_selection.max_threads,
                              "Thread limit for plan estimates; 0 uses the oneTBB default");
    methods->add_option("method", method_info, "Show details for a method ID");
    auto* parameter_set_option = parameters->add_option("parameter-set", parameter_set_info,
                                                        "Show details for a parameter-set ID");
    auto* parameter_method_option =
        parameters->add_option("--method", parameter_method, "Limit results to a method ID");
    parameter_set_option->excludes(parameter_method_option);
    parameters
        ->add_flag("--keys", parameter_keys, "List atom and bond parameter keys in matching order")
        ->needs(parameter_set_option);
    const auto argc = static_cast<int>(arguments.size());
    auto* const argv = arguments.data();
    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError& error) {
        return app.exit(error) == 0 ? 0 : 2;
    }

    if (*methods) {
        chargefw::cli::print_methods(method_info);
        return 0;
    }
    if (*parameters) {
        chargefw::cli::print_parameter_sets(parameter_set_info, parameter_method, parameter_keys);
        return 0;
    }
    if (*inspect) {
        chargefw::cli::print_inspection(chargefw::cli::import_input(inspect_input));
        return 0;
    }
    if (*applicability) {
        auto request = [&] {
            auto imported = chargefw::cli::import_input(applicability_input);
            return chargefw::cli::make_request(std::move(imported.molecules),
                                               applicability_selection, imported.records);
        }();
        chargefw::cli::print_applicability(chargefw::calculation::assess(std::move(request)),
                                           applicability_selection.max_threads);
        return 0;
    }
    auto run =
        chargefw::cli::CalculationRun{.metrics = {}, .started = std::chrono::steady_clock::now()};
    run.metrics.started_at = chargefw::cli::utc_timestamp();
    const auto parsing_started = std::chrono::steady_clock::now();
    auto imported = chargefw::cli::import_input(calculate_input);
    run.metrics.parsing_seconds =
        std::chrono::duration<double>{std::chrono::steady_clock::now() - parsing_started}.count();
    auto records = std::move(imported.records);
    auto request =
        chargefw::cli::make_request(std::move(imported.molecules), calculate_selection, records);
    const auto requested_provenance =
        chargefw::cli::make_requested_provenance(request, calculate_selection.max_threads);
    const auto max_threads = calculate_selection.max_threads;
    const auto progress_observer = TerminalProgressObserver{};
    const auto interrupt_handler = ScopedInterruptHandler{};
    const auto interruptible_observer =
        InterruptibleCalculationObserver{progress ? std::addressof(progress_observer) : nullptr};
    auto result = [&]() -> chargefw::calculation::ExecutionResult {
        try {
            auto assessment = chargefw::calculation::assess(std::move(request));
            const auto* plan = assessment.default_plan();
            if (plan == nullptr) {
                return chargefw::calculation::calculate(assessment, max_threads,
                                                        interruptible_observer);
            }
            for (const auto& warning : plan->warnings()) {
                std::println(std::cerr, "Warning: {}", warning.message);
            }
            return chargefw::calculation::calculate(assessment, *plan, max_threads,
                                                    interruptible_observer);
        } catch (const std::invalid_argument& error) {
            return {.status = chargefw::calculation::ExecutionStatus::invalid_input_or_request,
                    .charges = std::nullopt,
                    .failure_message = error.what()};
        } catch (const std::exception& error) {
            return {.status = chargefw::calculation::ExecutionStatus::numerical_failure,
                    .charges = std::nullopt,
                    .failure_message = error.what()};
        }
    }();
    run.metrics.applicability_seconds = result.metrics.applicability_seconds;
    run.metrics.computation_seconds = result.metrics.computation_seconds;
    return chargefw::cli::write_calculation_outputs(output_directory, calculate_input.path,
                                                    std::move(records), requested_provenance,
                                                    std::move(result), calculate_output, run);
}

} // namespace

auto main(int argc, char* argv[]) noexcept -> int {
    try {
        return run({argv, static_cast<std::size_t>(argc)});
    } catch (const std::invalid_argument& error) {
        std::print(std::cerr, "Invalid input or request: {}\n", error.what());
        return 2;
    } catch (const std::out_of_range& error) {
        std::print(std::cerr, "Invalid input or request: {}\n", error.what());
        return 2;
    } catch (const std::runtime_error& error) {
        std::print(std::cerr, "Invalid input or request: {}\n", error.what());
        return 2;
    } catch (const std::exception& error) {
        std::print(std::cerr, "Internal error: {}\n", error.what());
        return 1;
    }
}
