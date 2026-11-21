#include "support/test_exceptions.h"
#include "support/test_methods.h"

#include "calculation/target_execution.h"

#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/core/position.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method.h>
#include <chargefw/methods/method_metadata.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/methods/method_requirements.h>

#include <cstddef>
#include <exception>
#include <limits>
#include <optional>
#include <snitch/snitch.hpp>
#include <span>
#include <stdexcept>
#include <vector>

#include <oneapi/tbb/task_arena.h>

namespace calculation = chargefw::calculation;
namespace charges = chargefw::charges;
namespace core = chargefw::core;
namespace features = chargefw::features;
namespace methods = chargefw::methods;

namespace {

auto make_molecules(const std::size_t count) -> core::MoleculeCollection {
    auto molecules = std::vector<core::Molecule>{};
    molecules.reserve(count);
    for (std::size_t index = 0; index < count; ++index) {
        molecules.emplace_back(std::vector{core::Atom{6}}, std::vector<core::Bond>{},
                               std::vector<core::Conformer>{}, "target");
    }
    return core::MoleculeCollection{std::move(molecules)};
}

auto assert_worker_budget(const calculation::detail::ParallelizationLevel parallelization,
                          const std::size_t requested_threads, const std::size_t target_count,
                          const std::size_t expected_fragment_threads) -> void {
    const auto molecules = make_molecules(target_count);
    const features::PreparedMoleculeCollection prepared{molecules};
    const chargefw::test::StubMethod method{"target-execution"};
    const methods::ApplicableMethod selected{.method = &method, .parameter_set = nullptr};
    auto received_fragment_threads = std::vector<std::size_t>(target_count);
    const auto mode = parallelization == calculation::detail::ParallelizationLevel::targets
                          ? calculation::ExecutionMode::full
                          : calculation::ExecutionMode::cutoff;

    const auto result = calculation::detail::execute_calculation_targets(
        selected, prepared, mode, false, parallelization, requested_threads,
        calculation::default_calculation_observer(),
        [&](const features::PreparedMolecule& molecule,
            const chargefw::parameters::ParameterClassification*, const std::optional<std::size_t>,
            const std::size_t fragment_threads,
            const calculation::detail::ProgressContext& context) {
            received_fragment_threads[context.target_index] = fragment_threads;
            return charges::AtomicCharges{
                std::vector<double>(molecule.molecule().atom_count(), 0.0)};
        });

    REQUIRE(result.size() == target_count);
    for (std::size_t target_index = 0; target_index < target_count; ++target_index) {
        CHECK(result.assignment(target_index).target.molecule_index == target_index);
        CHECK(received_fragment_threads[target_index] == expected_fragment_threads);
    }
}

} // namespace

TEST_CASE("calculation targets distribute worker budget across parallelization levels",
          "[calculation][calculation-targets]") {
    for (const auto requested_threads : {std::size_t{0}, std::size_t{1}, std::size_t{2}}) {
        for (const auto target_count : {std::size_t{0}, std::size_t{1}, std::size_t{2}}) {
            assert_worker_budget(calculation::detail::ParallelizationLevel::targets,
                                 requested_threads, target_count, 1);
            assert_worker_budget(calculation::detail::ParallelizationLevel::fragments,
                                 requested_threads, target_count, requested_threads);
        }
    }
}

TEST_CASE("the default thread count stays within the allowed parallelism",
          "[calculation][calculation-targets]") {
    const auto allowed = calculation::detail::thread_count(std::numeric_limits<std::size_t>::max());
    const auto by_default = calculation::detail::thread_count(0);

    CHECK(by_default >= 1);
    CHECK(by_default <= allowed);
    CHECK(calculation::detail::thread_count(1) == 1);
}

TEST_CASE("an indexed loop runs its callable within the thread limit",
          "[calculation][calculation-targets]") {
    for (const auto requested_threads : {std::size_t{1}, std::size_t{2}}) {
        auto concurrency = std::vector<int>(3);
        calculation::detail::parallel_for_indexed(
            concurrency.size(), requested_threads, [&](const std::size_t index) {
                concurrency[index] = oneapi::tbb::this_task_arena::max_concurrency();
            });

        const auto expected =
            static_cast<int>(calculation::detail::thread_count(requested_threads));
        CHECK(concurrency == std::vector<int>(3, expected));
    }
}

TEST_CASE("calculation targets report method failures with their conformer",
          "[calculation][calculation-targets]") {
    const auto conformer = core::Conformer{{core::Position{}}};
    const auto molecules = core::MoleculeCollection{std::vector{
        core::Molecule{std::vector{core::Atom{6}}, {}, {conformer, conformer}, "pair"}}};
    const features::PreparedMoleculeCollection prepared{molecules};
    const chargefw::test::StubMethod method{"target-execution"};
    const methods::ApplicableMethod selected{.method = &method, .parameter_set = nullptr};
    // Fails on the second conformer only, with the given exception.
    const auto calculate = [&](const std::exception_ptr& error) {
        static_cast<void>(calculation::detail::execute_calculation_targets(
            selected, prepared, calculation::ExecutionMode::full, true,
            calculation::detail::ParallelizationLevel::targets, 1,
            calculation::default_calculation_observer(),
            [&](const features::PreparedMolecule& molecule,
                const chargefw::parameters::ParameterClassification*,
                const std::optional<std::size_t> conformer_index, const std::size_t,
                const calculation::detail::ProgressContext&) {
                if (conformer_index == 1) {
                    std::rethrow_exception(error);
                }
                return charges::AtomicCharges{
                    std::vector<double>(molecule.molecule().atom_count(), 0.0)};
            }));
    };

    const auto failed = chargefw::test::thrown<calculation::CalculationFailed>([&] {
        calculate(std::make_exception_ptr(
            methods::MethodFailure{methods::FailureKind::not_converged, "did not converge"}));
    });
    REQUIRE(failed.has_value());
    const auto* failure = &failed->failure();
    CHECK(failure->kind == methods::FailureKind::not_converged);
    CHECK(failure->message == "did not converge");
    CHECK(failure->molecule_index == 0);
    CHECK(failure->location == core::IssueLocation{.conformer_index = 1});

    // Anything else is an error in ChargeFW, described with its target. An invalid_argument is
    // itself a logic_error, so only the target shows that it was translated.
    CHECK_THROWS_AS(calculate(std::make_exception_ptr(std::runtime_error{"unexpected"})),
                    std::logic_error);
    CHECK_THROWS_MATCHES(calculate(std::make_exception_ptr(std::invalid_argument{"unexpected"})),
                         std::logic_error, snitch::matchers::with_what_contains{"conformer 2"});
}
