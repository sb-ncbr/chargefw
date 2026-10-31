#include "bindings.h"
#include "native_execution_result.h"
#include "native_input_metadata.h"
#include "native_parameter_catalog.h"

#include "adapters/fixed_ions.h"

#include <chargefw/adapters/charge_result.h>
#include <chargefw/calculation/assessment.h>
#include <chargefw/calculation/calculation.h>
#include <chargefw/calculation/execution_policy.h>
#include <chargefw/calculation/observer.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/methods/method.h>
#include <chargefw/methods/method_applicability.h>
#include <chargefw/methods/method_options.h>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include <cstddef>
#include <cstdint>
#include <map>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <variant>
#include <vector>

namespace nb = nanobind;

namespace chargefw::python {
namespace {

auto method_options(const nb::dict& values)
    -> std::unordered_map<std::string, methods::MethodOptions> {
    std::unordered_map<std::string, methods::MethodOptions> result;
    for (const auto& [method_value, options_value] : values) {
        const auto method_id = nb::cast<std::string>(method_value);
        const auto options = nb::cast<nb::dict>(options_value);
        std::unordered_map<std::string, methods::MethodOptionValue> native_values;
        for (const auto& [option_value, value] : options) {
            const auto option_id = nb::cast<std::string>(option_value);
            if (nb::isinstance<nb::bool_>(value)) {
                native_values.emplace(option_id, nb::cast<bool>(value));
            } else if (nb::isinstance<nb::int_>(value)) {
                native_values.emplace(option_id, nb::cast<int>(value));
            } else if (nb::isinstance<nb::float_>(value)) {
                native_values.emplace(option_id, nb::cast<double>(value));
            } else if (nb::isinstance<nb::str>(value)) {
                native_values.emplace(option_id, nb::cast<std::string>(value));
            } else {
                throw std::invalid_argument{
                    "method option values must be bool, int, float, or str"};
            }
        }
        result.emplace(method_id, methods::MethodOptions{std::move(native_values)});
    }
    return result;
}

[[nodiscard]] auto portable_id(const nb::handle value) -> adapters::PortableId {
    if (value.is_none()) {
        return {};
    }
    if (nb::isinstance<nb::str>(value)) {
        return adapters::PortableId{nb::cast<std::string>(value)};
    }
    if (nb::isinstance<nb::bool_>(value) || !nb::isinstance<nb::int_>(value)) {
        throw std::invalid_argument{"record ID must be a string, integer, or None"};
    }
    return adapters::PortableId{nb::cast<std::int64_t>(value)};
}

[[nodiscard]] auto portable_ids(const nb::handle value)
    -> std::optional<std::vector<adapters::PortableId>> {
    if (value.is_none()) {
        return std::nullopt;
    }
    auto result = std::vector<adapters::PortableId>{};
    const auto values = nb::cast<nb::sequence>(value);
    result.reserve(static_cast<std::size_t>(nb::len(values)));
    for (const auto& item : values) {
        result.push_back(portable_id(item));
    }
    return result;
}

auto prerequisite_issue(const methods::PrerequisiteIssue& issue) -> nb::dict {
    auto result = nb::dict{};
    result["kind"] = std::string{methods::to_string(issue.kind)};
    result["message"] = issue.message;
    result["molecule_index"] = issue.molecule_index;
    result["atom_index"] = issue.atom_index;
    result["bond_index"] = issue.bond_index;
    result["conformer_index"] = issue.conformer_index;
    return result;
}

auto execution_issue(const methods::ExecutionIssue& issue) -> nb::dict {
    auto result = nb::dict{};
    result["kind"] = std::string{methods::to_string(issue.kind)};
    result["message"] = issue.message;
    result["molecule_index"] = issue.molecule_index;
    return result;
}

auto execution_policy(const calculation::ExecutionPolicy& policy) -> nb::dict {
    auto result = nb::dict{};
    result["mode"] = std::string{calculation::to_string(policy.mode())};
    result["radius"] = policy.radius();
    return result;
}

auto method_option_value(const methods::MethodOptionValue& value) -> nb::object {
    return std::visit([](const auto& item) { return nb::cast(item); }, value);
}

auto method_options(const methods::MethodOptions& options) -> nb::dict {
    auto result = nb::dict{};
    for (const auto& [id, value] : options.values()) {
        result[id.c_str()] = method_option_value(value);
    }
    return result;
}

auto effective_calculation(const calculation::EffectiveCalculation& effective) -> nb::dict {
    auto result = nb::dict{};
    result["method_id"] = effective.method_id;
    result["parameter_set_id"] = effective.parameter_set_id;
    result["method_options"] = method_options(effective.method_options);
    result["execution_policy"] = execution_policy(effective.execution_policy);
    auto issues = nb::list{};
    for (const auto& issue : effective.execution_issues) {
        issues.append(execution_issue(issue));
    }
    result["execution_issues"] = std::move(issues);
    if (effective.fixed_ions.has_value()) {
        const auto& fixed_ions = *effective.fixed_ions;
        auto fixed_ions_value = nb::dict{};
        auto sources = nb::list{};
        for (const auto& source : fixed_ions.sources) {
            auto source_value = nb::dict{};
            source_value["molecule_index"] = source.molecule_index;
            source_value["atom_index"] = source.atom_index;
            source_value["charge"] = source.charge;
            sources.append(std::move(source_value));
        }
        fixed_ions_value["sources"] = std::move(sources);
        result["fixed_ions"] = std::move(fixed_ions_value);
    } else {
        result["fixed_ions"] = nb::none();
    }
    return result;
}

auto numpy_charge_values(const std::span<const double> values)
    -> nb::ndarray<nb::numpy, const double, nb::shape<-1>> {
    auto owned_values = std::make_unique<std::vector<double>>(values.begin(), values.end());
    auto owner = nb::capsule{owned_values.get(), [](void* pointer) noexcept {
                                 delete static_cast<std::vector<double>*>(pointer);
                             }};
    auto* released_values = owned_values.release();
    return {released_values->data(), {released_values->size()}, owner};
}

auto rejection(const calculation::Rejection& rejected) -> nb::dict {
    auto result = nb::dict{};
    result["method_id"] = rejected.method_id;
    result["parameter_set_id"] = rejected.parameter_set_id;
    result["execution_policy"] =
        rejected.policy.has_value() ? nb::cast(execution_policy(*rejected.policy)) : nb::none();
    auto prerequisite_issues = nb::list{};
    auto execution_issues = nb::list{};
    for (const auto& issue : rejected.issues) {
        if (const auto* prerequisite = std::get_if<methods::PrerequisiteIssue>(&issue)) {
            prerequisite_issues.append(prerequisite_issue(*prerequisite));
        } else {
            execution_issues.append(execution_issue(std::get<methods::ExecutionIssue>(issue)));
        }
    }
    result["prerequisite_issues"] = std::move(prerequisite_issues);
    result["execution_issues"] = std::move(execution_issues);
    return result;
}

auto execution_result(const calculation::ExecutionResult& value) -> nb::dict {
    auto result = nb::dict{};
    result["status"] = std::string{calculation::to_string(value.status)};
    auto rejections = nb::list{};
    for (const auto& rejected : value.rejections) {
        rejections.append(rejection(rejected));
    }
    result["rejections"] = std::move(rejections);
    result["failure_message"] = value.failure_message;
    auto metrics = nb::dict{};
    metrics["applicability_seconds"] = value.metrics.applicability_seconds;
    metrics["computation_seconds"] = value.metrics.computation_seconds;
    result["metrics"] = std::move(metrics);
    if (value.effective.has_value()) {
        result["effective"] = effective_calculation(*value.effective);
    } else {
        result["effective"] = nb::none();
    }

    if (!value.charges.has_value()) {
        result["charges"] = nb::none();
        return result;
    }

    const auto& charges = *value.charges;
    auto charge_set = nb::dict{};
    auto assignments = nb::list{};
    for (const auto& assignment : charges.assignments()) {
        auto item = nb::dict{};
        item["molecule_index"] = assignment.target.molecule_index;
        item["conformer_index"] = assignment.target.conformer_index;
        item["values"] = numpy_charge_values(assignment.charges.values());
        assignments.append(std::move(item));
    }
    charge_set["assignments"] = std::move(assignments);
    result["charges"] = std::move(charge_set);
    return result;
}

auto calculation_phase(const calculation::CalculationPhase phase) -> std::string_view {
    switch (phase) {
    case calculation::CalculationPhase::computation_started:
        return "computation_started";
    case calculation::CalculationPhase::computation_finished:
        return "computation_finished";
    case calculation::CalculationPhase::target_started:
        return "target_started";
    case calculation::CalculationPhase::target_finished:
        return "target_finished";
    case calculation::CalculationPhase::fragment_progress:
        return "fragment_progress";
    }
    throw std::invalid_argument{"unknown calculation phase"};
}

auto calculation_progress(const calculation::CalculationProgress& progress) -> nb::dict {
    auto result = nb::dict{};
    result["phase"] = std::string{calculation_phase(progress.phase)};
    result["mode"] = std::string{calculation::to_string(progress.mode)};
    result["method_id"] = std::string{progress.method_id};
    result["target_index"] = progress.target_index;
    result["target_count"] = progress.target_count;
    result["completed_fragment_count"] = progress.completed_fragment_count;
    result["fragment_count"] = progress.fragment_count;
    result["molecule_index"] = progress.molecule_index;
    result["conformer_index"] = progress.conformer_index;
    result["elapsed_seconds"] = progress.elapsed_seconds;
    return result;
}

class PythonCalculationObserver final : public calculation::CalculationObserver {
  public:
    explicit PythonCalculationObserver(nb::object observer) : observer_{std::move(observer)} {}

    void on_progress(const calculation::CalculationProgress& progress) const override {
        const auto acquire = nb::gil_scoped_acquire{};
        try {
            observer_.attr("on_progress")(calculation_progress(progress));
        } catch (nb::python_error& error) {
            error.discard_as_unraisable(observer_);
        }
    }

    [[nodiscard]] auto cancelled() const noexcept -> bool override {
        const auto acquire = nb::gil_scoped_acquire{};
        try {
            return nb::cast<bool>(observer_.attr("cancelled")());
        } catch (nb::python_error& error) {
            error.discard_as_unraisable(observer_);
        } catch (...) { // NOLINT(bugprone-empty-catch): cancellation callbacks cannot throw.
        }
        return false;
    }

  private:
    nb::object observer_;
};

class NativeAssessmentState {
  public:
    NativeAssessmentState(calculation::AssessmentResult assessment,
                          std::vector<adapters::ImportedMoleculeRecord> inputs,
                          adapters::RequestedCalculationProvenance requested,
                          const std::size_t max_threads)
        : assessment_{std::move(assessment)}, inputs_{std::move(inputs)},
          requested_{std::move(requested)}, max_threads_{max_threads} {}

    calculation::AssessmentResult assessment_;
    std::vector<adapters::ImportedMoleculeRecord> inputs_;
    adapters::RequestedCalculationProvenance requested_;
    std::size_t max_threads_ = 0;
};

auto execution_plan(const calculation::ExecutionPlan& plan) -> nb::dict {
    auto result = nb::dict{};
    result["method_id"] = std::string{plan.method().id()};
    result["parameter_set_id"] = plan.parameter_set() == nullptr
                                     ? nb::none()
                                     : nb::cast(std::string{plan.parameter_set()->id()});
    result["method_options"] = method_options(plan.method_options());
    result["execution_policy"] = execution_policy(plan.policy());
    auto warnings = nb::list{};
    for (const auto& warning : plan.warnings()) {
        warnings.append(execution_issue(warning));
    }
    result["warnings"] = std::move(warnings);
    return result;
}

class NativePlan {
  public:
    NativePlan(std::shared_ptr<NativeAssessmentState> state, const std::size_t index)
        : state_{std::move(state)}, index_{index} {}

    [[nodiscard]] auto report() const -> nb::dict {
        return execution_plan(state_->assessment_.plans()[index_]);
    }

    [[nodiscard]] auto calculate(const std::optional<std::size_t> max_threads,
                                 nb::object observer) const -> NativeExecutionResult {
        auto python_observer = std::optional<PythonCalculationObserver>{};
        const auto* native_observer = std::addressof(calculation::default_calculation_observer());
        if (!observer.is_none()) {
            python_observer.emplace(std::move(observer));
            native_observer = std::addressof(*python_observer);
        }
        auto result = [&] {
            nb::gil_scoped_release release;
            return calculation::calculate(state_->assessment_, state_->assessment_.plans()[index_],
                                          max_threads.value_or(state_->max_threads_),
                                          *native_observer);
        }();
        auto requested = state_->requested_;
        requested.max_threads = max_threads.value_or(state_->max_threads_);
        return NativeExecutionResult{adapters::make_charge_calculation_result(
            state_->inputs_, std::move(requested), std::move(result))};
    }

    [[nodiscard]] auto estimate(const std::optional<std::size_t> max_threads) const
        -> std::pair<double, double> {
        const auto estimate =
            calculation::estimate(state_->assessment_, state_->assessment_.plans()[index_],
                                  max_threads.value_or(state_->max_threads_));
        return {estimate.memory_bytes, estimate.operations};
    }

  private:
    std::shared_ptr<NativeAssessmentState> state_;
    std::size_t index_ = 0;
};

class NativeAssessment {
  public:
    NativeAssessment(calculation::AssessmentResult assessment,
                     std::vector<adapters::ImportedMoleculeRecord> inputs,
                     adapters::RequestedCalculationProvenance requested,
                     const std::size_t max_threads)
        : state_{std::make_shared<NativeAssessmentState>(std::move(assessment), std::move(inputs),
                                                         std::move(requested), max_threads)} {}

    [[nodiscard]] auto report() const -> nb::dict {
        auto result = nb::dict{};
        auto rejections = nb::list{};
        for (const auto& rejected : state_->assessment_.rejections()) {
            rejections.append(rejection(rejected));
        }
        result["rejections"] = std::move(rejections);
        result["applicability_seconds"] = state_->assessment_.applicability_seconds();
        return result;
    }

    [[nodiscard]] auto plans() const -> nb::list {
        auto result = nb::list{};
        for (std::size_t index = 0; index < state_->assessment_.plans().size(); ++index) {
            result.append(NativePlan{state_, index});
        }
        return result;
    }

    [[nodiscard]] auto calculate_default(nb::object observer) const -> NativeExecutionResult {
        auto python_observer = std::optional<PythonCalculationObserver>{};
        const auto* native_observer = std::addressof(calculation::default_calculation_observer());
        if (!observer.is_none()) {
            python_observer.emplace(std::move(observer));
            native_observer = std::addressof(*python_observer);
        }
        auto result = [&] {
            nb::gil_scoped_release release;
            return calculation::calculate(state_->assessment_, state_->max_threads_,
                                          *native_observer);
        }();
        return NativeExecutionResult{adapters::make_charge_calculation_result(
            state_->inputs_, state_->requested_, std::move(result))};
    }

  private:
    std::shared_ptr<NativeAssessmentState> state_;
};

auto make_assessment(const nb::sequence& molecules, const nb::sequence& input_metadata,
                     const nb::sequence& identities, const nb::sequence& caller_atom_ids,
                     std::string molecule_collection_name, const NativeParameterCatalog& catalog,
                     std::optional<std::string> method_id,
                     std::optional<std::string> parameter_set_id, const std::string& method_scope,
                     const nb::dict& options, const bool permissive_types,
                     const std::string& execution, const std::optional<double> radius,
                     const std::optional<std::size_t> cutoff_threshold,
                     const std::optional<std::size_t> cover_threshold,
                     const std::size_t max_threads, std::vector<std::string> fixed_ions)
    -> NativeAssessment {
    // Convert Python values to stable native references before releasing the GIL. The remaining
    // work copies native-owned input values and prepares the assessment without accessing Python
    // objects.
    const auto molecule_count = static_cast<std::size_t>(nb::len(molecules));
    if (static_cast<std::size_t>(nb::len(input_metadata)) != molecule_count ||
        static_cast<std::size_t>(nb::len(identities)) != molecule_count ||
        static_cast<std::size_t>(nb::len(caller_atom_ids)) != molecule_count) {
        throw std::invalid_argument{"input record metadata count does not match molecules"};
    }
    auto source_molecules = std::vector<const core::Molecule*>{};
    auto source_metadata = std::vector<const NativeInputMetadata*>{};
    auto source_identities = std::vector<adapters::MoleculeRecordIdentity>{};
    auto source_caller_atom_ids = std::vector<std::optional<std::vector<adapters::PortableId>>>{};
    source_molecules.reserve(molecule_count);
    source_metadata.reserve(molecule_count);
    source_identities.reserve(molecule_count);
    source_caller_atom_ids.reserve(molecule_count);
    for (std::size_t index = 0; index < molecule_count; ++index) {
        source_molecules.push_back(&nb::cast<const core::Molecule&>(molecules[index]));
        source_metadata.push_back(
            input_metadata[index].is_none()
                ? nullptr
                : &nb::cast<const NativeInputMetadata&>(input_metadata[index]));
        const auto identity = nb::cast<nb::sequence>(identities[index]);
        if (static_cast<std::size_t>(nb::len(identity)) != 3) {
            throw std::invalid_argument{"input record identity must contain three values"};
        }
        source_identities.push_back(
            adapters::MoleculeRecordIdentity{.source = nb::cast<std::string>(identity[0]),
                                             .record_index = nb::cast<std::size_t>(identity[1]),
                                             .record_id = portable_id(identity[2])});
        source_caller_atom_ids.push_back(portable_ids(caller_atom_ids[index]));
    }
    auto native_method_options = method_options(options);
    auto requested_options = std::map<std::string, methods::MethodOptions>{};
    for (const auto& [id, values] : native_method_options) {
        requested_options.emplace(id, values);
    }
    auto requested =
        adapters::RequestedCalculationProvenance{.method_id = method_id,
                                                 .parameter_set_id = parameter_set_id,
                                                 .method_scope = method_scope,
                                                 .permissive_types = permissive_types,
                                                 .cutoff_atom_threshold = cutoff_threshold,
                                                 .cover_atom_threshold = cover_threshold,
                                                 .max_threads = max_threads,
                                                 .execution_kind = execution,
                                                 .execution_radius = radius,
                                                 .method_options = std::move(requested_options)};
    nb::gil_scoped_release release;
    auto inputs = std::vector<adapters::ImportedMoleculeRecord>{};
    inputs.reserve(source_molecules.size());
    auto owned_molecules = std::vector<core::Molecule>{};
    owned_molecules.reserve(source_molecules.size());
    for (std::size_t index = 0; index < source_molecules.size(); ++index) {
        auto molecule = *source_molecules[index];
        inputs.push_back(
            source_metadata[index] == nullptr
                ? adapters::ImportedMoleculeRecord{.molecule = molecule,
                                                   .identity = std::move(source_identities[index]),
                                                   .caller_atom_ids =
                                                       std::move(source_caller_atom_ids[index]),
                                                   .diagnostics = {},
                                                   .import_metadata = std::nullopt}
                : source_metadata[index]->make_record(molecule));
        owned_molecules.push_back(std::move(molecule));
    }
    auto resolved_fixed_ions = adapters::detail::resolve_fixed_ions(inputs, fixed_ions);
    auto request = calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::move(owned_molecules),
                                              std::move(molecule_collection_name)},
        .parameter_sets = catalog.parameter_sets(),
        .method_id = std::move(method_id),
        .parameter_set_id = std::move(parameter_set_id),
        .method_scope = calculation::method_scope_from_string(method_scope),
        .method_options = std::move(native_method_options),
        .classification_options = {.permissive_types = permissive_types},
        .execution_selection =
            calculation::ExecutionSelection{
                calculation::execution_selection_kind_from_string(execution), radius},
        .resource_policy = {.cutoff_atom_threshold = cutoff_threshold,
                            .cover_atom_threshold = cover_threshold},
        .fixed_ions = resolved_fixed_ions.sources.empty()
                          ? std::nullopt
                          : std::optional{std::move(resolved_fixed_ions)},
    };
    return NativeAssessment{calculation::assess(std::move(request)), std::move(inputs),
                            std::move(requested), max_threads};
}

} // namespace

void bind_calculation(nb::module_& module) {
    module.def("_fixed_ion_names", [](const bool common_only) {
        nb::list names;
        for (const auto& name : adapters::detail::fixed_ion_names(common_only)) {
            names.append(nb::cast(name));
        }
        return nb::tuple{names};
    });
    nb::class_<NativeExecutionResult>(module, "_NativeExecutionResult")
        .def("report", [](const NativeExecutionResult& value) {
            return execution_result(value.result().execution());
        });
    nb::class_<NativePlan>(module, "_NativePlan")
        .def("report", &NativePlan::report)
        .def("calculate", &NativePlan::calculate, nb::arg("max_threads") = nb::none(),
             nb::arg("observer") = nb::none())
        .def("estimate", &NativePlan::estimate, nb::arg("max_threads") = nb::none());
    nb::class_<NativeAssessment>(module, "_NativeAssessment")
        .def("report", &NativeAssessment::report)
        .def("plans", &NativeAssessment::plans)
        .def("calculate_default", &NativeAssessment::calculate_default,
             nb::arg("observer") = nb::none());

    module.def("_make_assessment", &make_assessment, nb::arg("molecules"),
               nb::arg("input_metadata"), nb::arg("identities"), nb::arg("caller_atom_ids"),
               nb::arg("molecule_collection_name"), nb::arg("catalog"), nb::arg("method_id"),
               nb::arg("parameter_set_id"), nb::arg("method_scope"), nb::arg("method_options"),
               nb::arg("permissive_types"), nb::arg("execution"), nb::arg("radius"),
               nb::arg("cutoff_threshold"), nb::arg("cover_threshold"), nb::arg("max_threads"),
               nb::arg("fixed_ions") = std::vector<std::string>{});
}

} // namespace chargefw::python
