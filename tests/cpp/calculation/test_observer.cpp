#include "support/test_calculation.h"
#include "support/test_methods.h"
#include "support/test_molecules.h"
#include "support/test_parameters.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <chargefw/calculation/calculation.h>
#include <chargefw/calculation/observer.h>
#include <chargefw/core/atom.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method.h>
#include <chargefw/methods/method_applicability.h>
#include <chargefw/methods/method_metadata.h>
#include <chargefw/parameters/models/atom_parameters.h>
#include <chargefw/parameters/models/common_parameters.h>
#include <chargefw/parameters/models/parameter_set.h>
#include <chargefw/parameters/models/parameter_set_metadata.h>
#include <cstddef>
#include <cstdint>
#include <mutex>
#include <new>
#include <optional>
#include <snitch/snitch.hpp>
#include <stdexcept>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

namespace calculation = chargefw::calculation;
namespace charges = chargefw::charges;
namespace core = chargefw::core;
namespace features = chargefw::features;
namespace methods = chargefw::methods;

using chargefw::test::calculate_application;

namespace {

struct RecordedProgress {
    calculation::CalculationPhase phase{};
    calculation::ExecutionMode mode{};
    std::string method_id;
    std::size_t target_index{};
    std::size_t target_count{};
    std::size_t completed_fragment_count{};
    std::size_t fragment_count{};
    std::size_t molecule_index{};
    std::optional<std::size_t> conformer_index{};
    double elapsed_seconds{};
};

[[nodiscard]] auto snapshot(const calculation::CalculationProgress& progress) -> RecordedProgress {
    return RecordedProgress{.phase = progress.phase,
                            .mode = progress.mode,
                            .method_id = std::string{progress.method_id},
                            .target_index = progress.target_index,
                            .target_count = progress.target_count,
                            .completed_fragment_count = progress.completed_fragment_count,
                            .fragment_count = progress.fragment_count,
                            .molecule_index = progress.molecule_index,
                            .conformer_index = progress.conformer_index,
                            .elapsed_seconds = progress.elapsed_seconds};
}

// Thread-safe recording observer. Captures owned progress snapshots for later assertions.
class RecordingObserver : public calculation::CalculationObserver {
  public:
    void on_progress(const calculation::CalculationProgress& progress) const override {
        const std::scoped_lock lock{mutex_};
        events_.push_back(snapshot(progress));
    }

    [[nodiscard]] auto events() const -> std::vector<RecordedProgress> {
        const std::scoped_lock lock{mutex_};
        return events_;
    }

  private:
    mutable std::mutex mutex_;
    mutable std::vector<RecordedProgress> events_;
};

// Observer that records progress and requests cancellation after the first event in the given
// phase. Counting cancelled() calls after the request shows when a later loop iteration, rather
// than only the target boundary, observed the cancellation.
class CancellingObserver final : public calculation::CalculationObserver {
  public:
    explicit CancellingObserver(const calculation::CalculationPhase cancel_phase =
                                    calculation::CalculationPhase::target_started) noexcept
        : cancel_phase_{cancel_phase} {}

    void on_progress(const calculation::CalculationProgress& progress) const override {
        const std::scoped_lock lock{mutex_};
        events_.push_back(snapshot(progress));
        if (progress.phase == cancel_phase_) {
            cancel_.store(true, std::memory_order_relaxed);
        }
    }

    [[nodiscard]] auto cancelled() const noexcept -> bool override {
        if (cancel_.load(std::memory_order_relaxed)) {
            cancellation_checks_after_request_.fetch_add(1, std::memory_order_relaxed);
            return true;
        }
        return false;
    }

    [[nodiscard]] auto events() const -> std::vector<RecordedProgress> {
        const std::scoped_lock lock{mutex_};
        return events_;
    }

    [[nodiscard]] auto cancellation_checks_after_request() const noexcept -> std::size_t {
        return cancellation_checks_after_request_.load(std::memory_order_relaxed);
    }

  private:
    calculation::CalculationPhase cancel_phase_;
    mutable std::mutex mutex_;
    mutable std::vector<RecordedProgress> events_;
    mutable std::atomic<bool> cancel_{false};
    mutable std::atomic<std::size_t> cancellation_checks_after_request_{0};
};

[[nodiscard]] auto fragment_progress_events(const RecordingObserver& observer)
    -> std::vector<RecordedProgress> {
    auto result = std::vector<RecordedProgress>{};
    for (const auto& event : observer.events()) {
        if (event.phase == calculation::CalculationPhase::fragment_progress) {
            result.push_back(event);
        }
    }
    return result;
}

class ThrowOnEveryCallbackObserver final : public calculation::CalculationObserver {
  public:
    void on_progress(const calculation::CalculationProgress& /*progress*/) const override {
        callbacks_.fetch_add(1, std::memory_order_relaxed);
        throw std::runtime_error{"observer failure"};
    }

    [[nodiscard]] auto callbacks() const noexcept -> std::size_t {
        return callbacks_.load(std::memory_order_relaxed);
    }

  private:
    mutable std::atomic<std::size_t> callbacks_{0};
};

enum class DirectFailure : std::uint8_t { none, logic_error, allocation };

class DirectTestMethod final : public chargefw::test::StubMethod {
  public:
    explicit DirectTestMethod(const DirectFailure failure = DirectFailure::none)
        : StubMethod{"direct-test",
                     {.coordinates = true,
                      .resources = {.supports_cutoff = true,
                                    .supports_cover = true,
                                    .reduced_charge_policy =
                                        methods::ReducedChargePolicy::zero_components}}},
          failure_{failure} {}

    [[nodiscard]] auto calculate(const methods::CalculationInput& input) const
        -> charges::AtomicCharges override {
        if (failure_ == DirectFailure::logic_error) {
            throw std::logic_error{"direct observer test failure"};
        }
        if (failure_ == DirectFailure::allocation) {
            throw std::bad_alloc{};
        }
        return StubMethod::calculate(input);
    }

  private:
    DirectFailure failure_;
};

auto assert_single_terminal_fragment_progress(const std::vector<RecordedProgress>& events,
                                              const std::size_t fragment_count) -> void {
    REQUIRE(!events.empty());
    CHECK(std::count_if(events.begin(), events.end(), [fragment_count](const auto& event) {
              return event.completed_fragment_count == fragment_count &&
                     event.fragment_count == fragment_count;
          }) == 1);
    for (const auto& event : events) {
        CHECK(event.completed_fragment_count > 0);
        CHECK(event.completed_fragment_count <= event.fragment_count);
        CHECK(event.fragment_count == fragment_count);
    }
    for (std::size_t index = 1; index < events.size(); ++index) {
        CHECK(events[index - 1].completed_fragment_count < events[index].completed_fragment_count);
    }
}

auto make_separated_waters() -> core::Molecule {
    return core::Molecule{std::vector{core::Atom{8}, core::Atom{1}, core::Atom{1}, core::Atom{8},
                                      core::Atom{1}, core::Atom{1}},
                          std::vector{core::Bond{0, 1, core::BondOrder::SINGLE},
                                      core::Bond{0, 2, core::BondOrder::SINGLE},
                                      core::Bond{3, 4, core::BondOrder::SINGLE},
                                      core::Bond{3, 5, core::BondOrder::SINGLE}},
                          {core::Conformer{{{0.0, 0.0, 0.0},
                                            {0.96, 0.0, 0.0},
                                            {-0.24, 0.93, 0.0},
                                            {10.0, 0.0, 0.0},
                                            {10.96, 0.0, 0.0},
                                            {9.76, 0.93, 0.0}}}},
                          "separated-waters"};
}

auto make_many_separated_waters() -> core::Molecule {
    auto atoms = std::vector<core::Atom>{};
    auto bonds = std::vector<core::Bond>{};
    auto positions = std::vector<core::Position>{};
    constexpr std::size_t water_count = 12;
    atoms.reserve(water_count * 3);
    bonds.reserve(water_count * 2);
    positions.reserve(water_count * 3);

    for (std::size_t water_index = 0; water_index < water_count; ++water_index) {
        const auto atom_index = water_index * 3;
        const auto x = static_cast<double>(water_index) * 10.0;
        atoms.insert(atoms.end(), {core::Atom{8}, core::Atom{1}, core::Atom{1}});
        bonds.emplace_back(atom_index, atom_index + 1, core::BondOrder::SINGLE);
        bonds.emplace_back(atom_index, atom_index + 2, core::BondOrder::SINGLE);
        positions.insert(positions.end(),
                         {{x, 0.0, 0.0}, {x + 0.96, 0.0, 0.0}, {x - 0.24, 0.93, 0.0}});
    }

    return core::Molecule{std::move(atoms),
                          std::move(bonds),
                          {core::Conformer{std::move(positions)}},
                          "many-separated-waters"};
}

auto assert_computation_boundary(const std::vector<RecordedProgress>& events,
                                 const calculation::ExecutionMode mode) -> void {
    REQUIRE(!events.empty());
    CHECK(events.front().phase == calculation::CalculationPhase::computation_started);
    CHECK(events.back().phase == calculation::CalculationPhase::computation_finished);
    CHECK(events.front().mode == mode);
    CHECK(events.back().mode == mode);
    CHECK(std::count_if(events.begin(), events.end(), [](const auto& event) {
              return event.phase == calculation::CalculationPhase::computation_started;
          }) == 1);
    CHECK(std::count_if(events.begin(), events.end(), [](const auto& event) {
              return event.phase == calculation::CalculationPhase::computation_finished;
          }) == 1);
}

constexpr std::array all_modes{calculation::ExecutionMode::full, calculation::ExecutionMode::cutoff,
                               calculation::ExecutionMode::cover};
constexpr std::array reduced_modes{calculation::ExecutionMode::cutoff,
                                   calculation::ExecutionMode::cover};

// Explicit selection of one execution mode; reduced modes use the minimum radius.
[[nodiscard]] auto select_mode(const calculation::ExecutionMode mode)
    -> calculation::ExecutionSelection {
    switch (mode) {
    case calculation::ExecutionMode::full:
        return calculation::ExecutionSelection{calculation::ExecutionSelectionKind::full};
    case calculation::ExecutionMode::cutoff:
        return calculation::ExecutionSelection{calculation::ExecutionSelectionKind::cutoff,
                                               calculation::minimum_reduced_radius};
    case calculation::ExecutionMode::cover:
        return calculation::ExecutionSelection{calculation::ExecutionSelectionKind::cover,
                                               calculation::minimum_reduced_radius};
    }
    throw std::logic_error{"unknown execution mode"};
}

[[nodiscard]] auto make_direct_policy(const calculation::ExecutionMode mode)
    -> calculation::ExecutionPolicy {
    return mode == calculation::ExecutionMode::full
               ? calculation::ExecutionPolicy{}
               : calculation::ExecutionPolicy{mode, calculation::minimum_reduced_radius};
}

// QEq parameters whose zero hydrogen hardness makes the charge solve fail numerically.
auto make_invalid_qeq_parameters() -> chargefw::parameters::ParameterSet {
    return chargefw::parameters::ParameterSet{
        chargefw::parameters::ParameterSetMetadata{
            .id = "invalid-qeq", .method_id = "qeq", .name = "Invalid QEq"},
        {},
        chargefw::parameters::AtomParameters{
            {{.key = chargefw::test::plain_atom_key(1),
              .parameters = {{.name = "electronegativity", .value = 4.5280},
                             {.name = "hardness", .value = 0.0}}},
             {.key = chargefw::test::plain_atom_key(8),
              .parameters = {{.name = "electronegativity", .value = 8.741},
                             {.name = "hardness", .value = 13.364}}}}}};
}

auto make_formal_water_request() -> calculation::AssessmentRequest {
    return calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .method_id = "formal",
        .execution_selection =
            calculation::ExecutionSelection{calculation::ExecutionSelectionKind::full}};
}

} // namespace

TEST_CASE("observer emits ordered computation and target phases", "[calculation][observer]") {
    const auto observer = RecordingObserver{};
    const auto result = calculate_application(make_formal_water_request(), observer);

    REQUIRE(result.calculated());
    CHECK(!result.cancelled());

    const auto events = observer.events();
    assert_computation_boundary(events, calculation::ExecutionMode::full);
    CHECK(events.front().method_id == std::string_view{"formal"});
    CHECK(events.back().method_id == std::string_view{"formal"});

    // Target phases occur between the computation boundaries.
    const auto computation_end = events.end() - 1;
    const auto target_started = std::find_if(events.begin(), computation_end, [](const auto& e) {
        return e.phase == calculation::CalculationPhase::target_started;
    });
    CHECK(target_started != computation_end);
    const auto target_finished = std::find_if(target_started, computation_end, [](const auto& e) {
        return e.phase == calculation::CalculationPhase::target_finished;
    });
    CHECK(target_finished != computation_end);
}

TEST_CASE("cancellation produces a terminal observer event", "[calculation][observer]") {
    const auto observer = CancellingObserver{};
    const auto result = calculate_application(make_formal_water_request(), observer);

    CHECK(result.cancelled());
    CHECK(!result.calculated());
    CHECK(result.status == calculation::ExecutionStatus::cancelled);
    CHECK(!result.charges.has_value());
    REQUIRE(result.effective.has_value());
    CHECK(result.effective->method_id == "formal");
    assert_computation_boundary(observer.events(), calculation::ExecutionMode::full);
}

TEST_CASE("fixed-charge cancellation retains provenance and permits plan reuse",
          "[calculation][observer]") {
    const auto molecule =
        core::Molecule{std::vector{core::Atom{1, 0}, core::Atom{12, 2}, core::Atom{8, 0}},
                       {},
                       {core::Conformer{{{0.0, 0.0, 0.0}, {0.0, 3.0, 0.0}, {5.0, 0.0, 0.0}}}},
                       "cancelled-fixed-charge-eem"};
    for (const auto mode : all_modes) {
        CAPTURE(calculation::to_string(mode));
        const auto observer =
            CancellingObserver{mode == calculation::ExecutionMode::full
                                   ? calculation::CalculationPhase::target_started
                                   : calculation::CalculationPhase::fragment_progress};
        const auto assessment = calculation::assess(calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::vector{molecule}},
            .parameter_sets = {chargefw::test::make_eem_ho_parameters(2.0)},
            .method_id = "eem",
            .execution_selection = select_mode(mode),
            .fixed_ions = calculation::FixedIons{
                .sources = {{.molecule_index = 0, .atom_index = 1, .charge = 0.4}}}});
        REQUIRE(assessment.plans().size() == 1);
        CHECK(assessment.plans()[0].policy().mode() == mode);
        const auto result = calculation::calculate(assessment, 1, observer);

        CHECK(result.status == calculation::ExecutionStatus::cancelled);
        CHECK_FALSE(result.calculated());
        CHECK_FALSE(result.charges.has_value());
        REQUIRE(result.effective.has_value());
        REQUIRE(result.effective->fixed_ions.has_value());
        const auto& provenance = *result.effective->fixed_ions;
        CHECK(provenance.sources.size() == 1);
        CHECK(provenance.sources[0].charge == 0.4);

        const auto events = observer.events();
        assert_computation_boundary(events, mode);
        CHECK(std::count_if(events.begin(), events.end(), [](const auto& event) {
                  return event.phase == calculation::CalculationPhase::target_started;
              }) == 1);
        if (mode == calculation::ExecutionMode::full) {
            CHECK(std::ranges::none_of(events, [](const auto& event) {
                return event.phase == calculation::CalculationPhase::fragment_progress;
            }));
        } else {
            const auto fragment_progress = std::ranges::find_if(events, [](const auto& event) {
                return event.phase == calculation::CalculationPhase::fragment_progress;
            });
            REQUIRE(fragment_progress != events.end());
            CHECK(fragment_progress->completed_fragment_count == 1);
            CHECK(fragment_progress->fragment_count == 2);
        }

        const auto repeated = calculation::calculate(assessment);
        REQUIRE(repeated.calculated());
        REQUIRE(repeated.effective->fixed_ions.has_value());
        CHECK(repeated.effective->fixed_ions->sources[0].charge == 0.4);
    }
}

TEST_CASE("solver failures finish observation with a numerical result in every mode",
          "[calculation][observer]") {
    for (const auto mode : all_modes) {
        CAPTURE(calculation::to_string(mode));
        const auto observer = RecordingObserver{};
        const auto result = calculate_application(
            calculation::AssessmentRequest{
                .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
                .parameter_sets = {make_invalid_qeq_parameters()},
                .method_id = "qeq",
                .parameter_set_id = "invalid-qeq",
                .execution_selection = select_mode(mode)},
            observer);

        CHECK(result.status == calculation::ExecutionStatus::numerical_failure);
        CHECK_FALSE(result.calculated());
        REQUIRE(result.failure_message.has_value());
        CHECK(result.failure_message->contains("method 'qeq', molecule 1 ('water'), conformer 1"));
        if (mode != calculation::ExecutionMode::full) {
            CHECK(result.failure_message->contains(std::string{calculation::to_string(mode)} +
                                                   " fragment around source atom 1 failed"));
        }

        const auto events = observer.events();
        assert_computation_boundary(events, mode);
        CHECK(events.back().method_id == std::string_view{"qeq"});
    }
}

TEST_CASE("serial cutoff execution emits aggregate fragment progress", "[calculation][observer]") {
    const auto observer = RecordingObserver{};
    const auto result = calculate_application(
        calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
            .execution_selection = select_mode(calculation::ExecutionMode::cutoff)},
        observer);

    REQUIRE(result.calculated());
    CHECK(!result.cancelled());

    const auto fragment_events = fragment_progress_events(observer);
    REQUIRE(!fragment_events.empty());
    assert_single_terminal_fragment_progress(fragment_events,
                                             fragment_events.front().fragment_count);
}

TEST_CASE("serial cover execution emits aggregate multi-pivot progress",
          "[calculation][observer]") {
    const auto observer = RecordingObserver{};
    const auto result = calculate_application(
        calculation::AssessmentRequest{
            .molecules = core::MoleculeCollection{std::vector{make_separated_waters()}},
            .execution_selection = select_mode(calculation::ExecutionMode::cover)},
        observer);

    REQUIRE(result.calculated());
    CHECK(!result.cancelled());

    const auto fragment_events = fragment_progress_events(observer);
    REQUIRE(!fragment_events.empty());
    assert_single_terminal_fragment_progress(fragment_events,
                                             fragment_events.front().fragment_count);
}

TEST_CASE("empty reduced targets do not emit fragment progress", "[calculation][observer]") {
    for (const auto mode : reduced_modes) {
        const auto observer = RecordingObserver{};
        const auto empty_molecule = core::Molecule{{}, {}, {core::Conformer{{}}}, "empty"};
        const auto result = calculate_application(
            calculation::AssessmentRequest{
                .molecules = core::MoleculeCollection{std::vector{empty_molecule}},
                .execution_selection = select_mode(mode)},
            observer);
        REQUIRE(result.calculated());
        CHECK(fragment_progress_events(observer).empty());
    }
}

TEST_CASE("parallel reduced targets each emit one terminal progress snapshot",
          "[calculation][observer]") {
    for (const auto mode : reduced_modes) {
        const auto observer = RecordingObserver{};
        const auto result = calculate_application(
            calculation::AssessmentRequest{
                .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water(),
                                                                  chargefw::test::make_water()}},
                .execution_selection = select_mode(mode)},
            observer, 2);
        CHECK(result.calculated());

        for (const auto target_index : {std::size_t{0}, std::size_t{1}}) {
            auto target_events = std::vector<RecordedProgress>{};
            for (const auto& event : fragment_progress_events(observer)) {
                if (event.target_index == target_index) {
                    target_events.push_back(event);
                }
            }
            REQUIRE(!target_events.empty());
            assert_single_terminal_fragment_progress(target_events,
                                                     target_events.front().fragment_count);
        }
    }
}

TEST_CASE("results and target events preserve source target identity in every execution mode",
          "[calculation][observer]") {
    const auto expected_targets = std::vector<charges::ChargeTarget>{
        {.molecule_index = 0, .conformer_index = 0},
        {.molecule_index = 0, .conformer_index = 1},
        {.molecule_index = 1, .conformer_index = 0},
    };
    for (const auto mode : all_modes) {
        for (const auto max_threads : {std::size_t{1}, std::size_t{2}}) {
            CAPTURE(calculation::to_string(mode), max_threads);
            const auto observer = RecordingObserver{};
            const auto result = calculate_application(
                calculation::AssessmentRequest{
                    .molecules = core::MoleculeCollection{std::vector{
                        chargefw::test::make_two_conformer_water(), chargefw::test::make_water()}},
                    .method_id = "eqeq",
                    .execution_selection = select_mode(mode)},
                observer, max_threads);
            REQUIRE(result.calculated());
            REQUIRE(result.charges->size() == expected_targets.size());
            for (std::size_t index = 0; index < expected_targets.size(); ++index) {
                CHECK(result.charges->assignment(index).target.molecule_index ==
                      expected_targets[index].molecule_index);
                CHECK(result.charges->assignment(index).target.conformer_index ==
                      expected_targets[index].conformer_index);
            }

            // Parallel targets may start in any order, but each event carries its own identity.
            auto started = std::vector<bool>(expected_targets.size(), false);
            for (const auto& event : observer.events()) {
                if (event.phase != calculation::CalculationPhase::target_started) {
                    continue;
                }
                REQUIRE(event.target_index < expected_targets.size());
                CHECK_FALSE(started[event.target_index]);
                started[event.target_index] = true;
                CHECK(event.target_count == expected_targets.size());
                CHECK(event.molecule_index == expected_targets[event.target_index].molecule_index);
                CHECK(event.conformer_index ==
                      expected_targets[event.target_index].conformer_index);
            }
            CHECK(std::ranges::all_of(started, [](const bool value) { return value; }));
        }
    }
}

TEST_CASE("unsupported explicit execution is rejected before calculation observation begins",
          "[calculation][observer]") {
    const auto observer = RecordingObserver{};
    const auto assessment = calculation::assess(calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}},
        .method_id = "formal",
        .execution_selection = select_mode(calculation::ExecutionMode::cover)});

    CHECK(assessment.plans().empty());
    REQUIRE(assessment.rejections().size() == 1);
    REQUIRE(assessment.rejections()[0].policy.has_value());
    CHECK(assessment.rejections()[0].policy->mode() == calculation::ExecutionMode::cover);
    REQUIRE(assessment.rejections()[0].issues.size() == 1);
    CHECK(std::get<methods::ExecutionIssue>(assessment.rejections()[0].issues[0]).kind ==
          methods::ExecutionIssueKind::unsupported_execution_mode);

    const auto result = calculation::calculate(assessment, 1, observer);
    CHECK(result.status == calculation::ExecutionStatus::no_executable_plan);
    CHECK_FALSE(result.calculated());
    CHECK(observer.events().empty());
}

TEST_CASE("observer callback failures do not alter calculation control flow",
          "[calculation][observer]") {
    const auto observer = ThrowOnEveryCallbackObserver{};
    const auto result = calculate_application(make_formal_water_request(), observer);
    REQUIRE(result.calculated());
    CHECK(!result.cancelled());
    CHECK(observer.callbacks() >= 4);
}

TEST_CASE("reduced execution observes cancellation after fragment progress",
          "[calculation][observer]") {
    for (const auto mode : reduced_modes) {
        for (const auto max_threads : {std::size_t{1}, std::size_t{2}}) {
            CAPTURE(calculation::to_string(mode), max_threads);
            const auto observer =
                CancellingObserver{calculation::CalculationPhase::fragment_progress};
            const auto result = calculate_application(
                calculation::AssessmentRequest{.molecules = core::MoleculeCollection{std::vector{
                                                   make_many_separated_waters()}},
                                               .execution_selection = select_mode(mode)},
                observer, max_threads);

            CHECK(result.cancelled());
            CHECK(!result.calculated());
            CHECK(result.status == calculation::ExecutionStatus::cancelled);
            CHECK(observer.cancellation_checks_after_request() > 0);

            const auto events = observer.events();
            assert_computation_boundary(events, mode);
            CHECK(std::ranges::any_of(events, [](const auto& event) {
                return event.phase == calculation::CalculationPhase::fragment_progress;
            }));
            CHECK(std::ranges::none_of(events, [](const auto& event) {
                return event.phase == calculation::CalculationPhase::target_finished;
            }));
        }
    }
}

TEST_CASE("direct calculation finishes observation after success, failure, and cancellation",
          "[calculation][observer]") {
    const auto molecules = core::MoleculeCollection{std::vector{chargefw::test::make_water()}};
    const auto prepared = features::PreparedMoleculeCollection{molecules};
    const auto method = DirectTestMethod{};
    const auto failing_method = DirectTestMethod{DirectFailure::logic_error};
    const auto allocating_method = DirectTestMethod{DirectFailure::allocation};

    for (const auto mode : all_modes) {
        CAPTURE(calculation::to_string(mode));
        const auto calculate_direct = [&](const methods::Method& selected_method,
                                          const calculation::CalculationObserver& observer,
                                          const std::size_t max_threads = 1) {
            const auto selected = methods::ApplicableMethod{.method = &selected_method};
            return calculation::calculate({.molecules = prepared,
                                           .selected = selected,
                                           .execution_policy = make_direct_policy(mode),
                                           .max_threads = max_threads,
                                           .observer = observer});
        };

        const auto observer = RecordingObserver{};
        CHECK(calculate_direct(method, observer).charges.size() == 1);
        assert_computation_boundary(observer.events(), mode);

        const auto failure_observer = RecordingObserver{};
        // The direct overload throws failures; it keeps their kind while adding the target or
        // fragment. The assessment-based overload reports this one as numerical_failure.
        CHECK_THROWS_AS(static_cast<void>(calculate_direct(failing_method, failure_observer)),
                        std::logic_error);
        assert_computation_boundary(failure_observer.events(), mode);

        // An allocation failure is not wrapped by the target or fragment handlers, also when
        // fragments or targets run in parallel, so callers can report it as running out of memory.
        for (const auto max_threads : {std::size_t{1}, std::size_t{4}}) {
            CAPTURE(max_threads);
            const auto allocation_observer = RecordingObserver{};
            CHECK_THROWS_AS(static_cast<void>(calculate_direct(allocating_method,
                                                               allocation_observer, max_threads)),
                            std::bad_alloc);
            assert_computation_boundary(allocation_observer.events(), mode);
        }

        const auto cancelling_observer = CancellingObserver{};
        CHECK_THROWS_AS(static_cast<void>(calculate_direct(method, cancelling_observer)),
                        calculation::CalculationCancelled);
        assert_computation_boundary(cancelling_observer.events(), mode);
    }
}
