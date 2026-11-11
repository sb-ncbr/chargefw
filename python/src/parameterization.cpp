#include "bindings.h"
#include "calculation.h"
#include "native_parameter_catalog.h"
#include "numpy_values.h"

#include <chargefw/calculation/assessment.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method.h>
#include <chargefw/methods/method_registry.h>
#include <chargefw/parameterization/normal_equations.h>
#include <chargefw/parameterization/parameter_entries.h>
#include <chargefw/parameterization/parameter_template.h>
#include <chargefw/parameterization/parameter_vector.h>
#include <chargefw/parameters/classification/parameter_classifier.h>
#include <chargefw/parameters/models/parameter_set.h>
#include <chargefw/parameters/models/parameter_set_metadata.h>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

namespace nb = nanobind;

namespace chargefw::python {
namespace {

// The slots of a parameter set, as (scope, entry index, name).
auto layout_list(const parameters::ParameterSet& parameter_set) -> nb::list {
    auto result = nb::list{};
    for (const auto& slot : parameterization::parameter_layout(parameter_set)) {
        result.append(nb::make_tuple(std::string{parameterization::to_string(slot.scope)},
                                     slot.entry_index, slot.name));
    }
    return result;
}

// The fitting view of one plan: its parameter values, the charges and normal equations at other
// values, and what the molecules show about its parameter entries. Shares the plan's assessment
// state, which it keeps alive.
class NativeEvaluator {
  public:
    explicit NativeEvaluator(const NativePlan& plan)
        : state_{plan.state()}, plan_{std::addressof(plan.plan())} {
        static_cast<void>(parameter_set());
    }

    // The slots of the plan's parameter set, as (scope, entry index, name).
    [[nodiscard]] auto layout() const -> nb::list {
        return layout_list(parameter_set());
    }

    [[nodiscard]] auto values() const -> nb::ndarray<nb::numpy, const double, nb::shape<-1>> {
        return numpy_values(parameterization::parameter_values(parameter_set()));
    }

    // Charges of every assignment, concatenated in assignment order, for candidate values in
    // layout order. Only the status, the charges with the size and molecule of each assignment, and
    // the failure message are returned: a fit evaluates many candidates and needs no other result
    // data.
    [[nodiscard]] auto calculate_values(const value_array& values,
                                        const std::optional<std::size_t> max_threads) const
        -> nb::dict {
        auto candidate_values = std::vector<double>(values.data(), values.data() + values.size());
        auto result = [&] {
            nb::gil_scoped_release release;
            return calculation::calculate(plan(), candidate(candidate_values),
                                          max_threads.value_or(state_->max_threads_));
        }();

        auto payload = nb::dict{};
        payload["status"] = std::string{calculation::to_string(result.status)};
        payload["failure_message"] =
            result.failure_message.has_value() ? nb::cast(*result.failure_message) : nb::none();
        if (!result.charges.has_value()) {
            payload["charges"] = nb::none();
            payload["assignment_sizes"] = nb::none();
            payload["assignment_molecules"] = nb::none();
            return payload;
        }
        auto charges = std::vector<double>{};
        auto assignment_sizes = nb::list{};
        auto assignment_molecules = nb::list{};
        for (const auto& assignment : result.charges->assignments()) {
            const auto assignment_values = assignment.charges.values();
            charges.insert(charges.end(), assignment_values.begin(), assignment_values.end());
            assignment_sizes.append(assignment_values.size());
            assignment_molecules.append(assignment.target.molecule_index);
        }
        payload["charges"] = numpy_values(charges);
        payload["assignment_sizes"] = std::move(assignment_sizes);
        payload["assignment_molecules"] = std::move(assignment_molecules);
        return payload;
    }

    // The normal equations of a least-squares fit at candidate values in layout order
    // (parameterization::normal_equations); "gauss_newton" is flattened row by row. A failing
    // calculation is reported in "status" and "failure_message", as by calculate_values().
    [[nodiscard]] auto normal_equations(const value_array& values, const value_array& references,
                                        const value_array& weights,
                                        const std::optional<std::size_t> max_threads) const
        -> nb::dict {
        auto candidate_values = std::vector<double>(values.data(), values.data() + values.size());
        auto result = std::optional<parameterization::NormalEquations>{};
        auto failure_message = std::optional<std::string>{};
        {
            nb::gil_scoped_release release;
            try {
                result = parameterization::normal_equations(
                    plan(), candidate(candidate_values),
                    std::span{references.data(), references.size()},
                    std::span{weights.data(), weights.size()},
                    max_threads.value_or(state_->max_threads_));
            } catch (const std::invalid_argument&) {
                throw;
            } catch (const std::logic_error& error) {
                failure_message = error.what();
            } catch (const std::runtime_error& error) {
                failure_message = error.what();
            }
        }

        auto payload = nb::dict{};
        payload["status"] = result.has_value() ? "success" : "numerical_failure";
        payload["failure_message"] =
            failure_message.has_value() ? nb::cast(*failure_message) : nb::none();
        if (!result.has_value()) {
            payload["charges"] = nb::none();
            payload["sum_of_squares"] = nb::none();
            payload["gauss_newton"] = nb::none();
            payload["gradient"] = nb::none();
            return payload;
        }
        payload["charges"] = numpy_values(result->charges);
        payload["sum_of_squares"] = result->sum_of_squares;
        payload["gauss_newton"] = numpy_values(result->gauss_newton);
        payload["gradient"] = numpy_values(result->gradient);
        return payload;
    }

    [[nodiscard]] auto exact_derivatives() const -> bool {
        return parameterization::has_exact_derivatives(plan().method());
    }

    // The number of calculated atoms matched to each atom entry and of bonds matched to each bond
    // entry of the plan's parameter set (parameterization::entry_counts()).
    [[nodiscard]] auto entry_counts() const
        -> std::pair<std::vector<std::size_t>, std::vector<std::size_t>> {
        auto counts = parameterization::entry_counts(plan());
        return {std::move(counts.atoms), std::move(counts.bonds)};
    }

    // Per molecule of the plan, the atom entry of every atom, -1 where none matches, as
    // parameterization::atom_entries() reports them.
    [[nodiscard]] auto atom_entries() const -> nb::list {
        auto result = nb::list{};
        for (const auto& molecule : parameterization::atom_entries(plan())) {
            auto entries = std::vector<std::int64_t>{};
            entries.reserve(molecule.size());
            for (const auto& entry : molecule) {
                entries.push_back(entry.has_value() ? static_cast<std::int64_t>(*entry) : -1);
            }
            result.append(numpy_integers(std::move(entries)));
        }
        return result;
    }

    // Per atom entry: (single_bonded, exclusive_bonds, bond_entries), as
    // parameterization::single_bond_entries() reports them.
    [[nodiscard]] auto single_bond_entries() const
        -> std::vector<std::tuple<bool, bool, std::vector<std::size_t>>> {
        auto result = std::vector<std::tuple<bool, bool, std::vector<std::size_t>>>{};
        for (auto& entry : parameterization::single_bond_entries(plan())) {
            result.emplace_back(entry.single_bonded, entry.exclusive_bonds,
                                std::move(entry.bond_entries));
        }
        return result;
    }

  private:
    [[nodiscard]] auto plan() const -> const calculation::ExecutionPlan& {
        return *plan_;
    }

    // The plan's parameter set with `values`, under an ID of its own.
    [[nodiscard]] auto candidate(const std::span<const double> values) const
        -> parameters::ParameterSet {
        const auto& assessed = parameter_set();
        auto metadata = parameters::ParameterSetMetadata{};
        metadata.id = std::string{assessed.id()} + ":candidate";
        return parameterization::with_parameter_values(assessed, values, std::move(metadata));
    }

    [[nodiscard]] auto parameter_set() const -> const parameters::ParameterSet& {
        const auto* parameter_set = plan().parameter_set();
        if (parameter_set == nullptr) {
            throw std::invalid_argument{"method '" + std::string{plan().method().id()} +
                                        "' has no parameters"};
        }
        return *parameter_set;
    }

    std::shared_ptr<NativeAssessmentState> state_;
    // Owned by the assessment in state_.
    const calculation::ExecutionPlan* plan_ = nullptr;
};

auto parameter_names(const std::vector<std::string_view>& names) -> std::vector<std::string> {
    return {names.begin(), names.end()};
}

// A catalog with a template for `method_id` built from the molecules: one zero-valued entry for
// every atom and bond type that occurs, with the parameter names the method requires.
auto make_parameter_template_catalog(const nb::sequence& molecules, const std::string& method_id,
                                     const std::string& atom_types, const std::string& bond_types,
                                     std::string id) -> NativeParameterCatalog {
    auto copies = std::vector<core::Molecule>{};
    copies.reserve(static_cast<std::size_t>(nb::len(molecules)));
    for (const auto molecule : molecules) {
        copies.push_back(nb::cast<const core::Molecule&>(molecule));
    }

    nb::gil_scoped_release release;
    const auto* method = methods::method_registry().find(method_id);
    if (method == nullptr) {
        throw std::invalid_argument{"method '" + method_id + "' is not registered"};
    }
    const auto requirements = method->requirements();
    if (!requirements.requires_parameters()) {
        throw std::invalid_argument{"method '" + method_id + "' has no parameters"};
    }
    const auto collection = core::MoleculeCollection{std::move(copies)};
    const auto prepared = features::PreparedMoleculeCollection{collection};
    auto metadata = parameters::ParameterSetMetadata{};
    metadata.id = std::move(id);
    metadata.method_id = method_id;
    auto parameter_sets = std::vector<parameters::ParameterSet>{};
    parameter_sets.push_back(parameterization::make_parameter_template(
        prepared,
        {.common = parameter_names(requirements.common_parameters),
         .atom = parameter_names(requirements.atom_parameters),
         .bond = parameter_names(requirements.bond_parameters)},
        parameters::atom_classification_kind_from_string(atom_types),
        parameters::bond_classification_kind_from_string(bond_types), std::move(metadata)));
    return NativeParameterCatalog{std::move(parameter_sets)};
}

// The slots and values of the parameter set `id`, in layout order.
auto parameter_set_values(const NativeParameterCatalog& catalog, const std::string& id)
    -> nb::tuple {
    const auto& parameter_set = catalog.find(id);
    return nb::make_tuple(layout_list(parameter_set),
                          numpy_values(parameterization::parameter_values(parameter_set)));
}

// A catalog with a copy of the parameter set `id`, with `values` in layout order and new metadata;
// the method stays the same.
auto parameter_set_with_values(const NativeParameterCatalog& catalog, const std::string& id,
                               const value_array& values, std::string new_id, std::string name,
                               std::string publication, std::string notes)
    -> NativeParameterCatalog {
    auto metadata = parameters::ParameterSetMetadata{};
    metadata.id = std::move(new_id);
    metadata.name = std::move(name);
    metadata.publication = std::move(publication);
    metadata.notes = std::move(notes);
    auto parameter_sets = std::vector<parameters::ParameterSet>{};
    parameter_sets.push_back(parameterization::with_parameter_values(
        catalog.find(id), std::span{values.data(), values.size()}, std::move(metadata)));
    return NativeParameterCatalog{std::move(parameter_sets)};
}

} // namespace

void bind_parameterization(nb::module_& module) {
    nb::class_<NativeEvaluator>(module, "_NativeEvaluator")
        .def(nb::init<const NativePlan&>(), nb::arg("plan"))
        .def("_layout", &NativeEvaluator::layout)
        .def("_values", &NativeEvaluator::values)
        .def("_calculate_values", &NativeEvaluator::calculate_values, nb::arg("values"),
             nb::arg("max_threads") = nb::none())
        .def("_entry_counts", &NativeEvaluator::entry_counts)
        .def("_atom_entries", &NativeEvaluator::atom_entries)
        .def("_single_bond_entries", &NativeEvaluator::single_bond_entries)
        .def("_normal_equations", &NativeEvaluator::normal_equations, nb::arg("values"),
             nb::arg("references"), nb::arg("weights"), nb::arg("max_threads") = nb::none())
        .def("_exact_derivatives", &NativeEvaluator::exact_derivatives);
    module.def("_make_parameter_template_catalog", &make_parameter_template_catalog,
               nb::arg("molecules"), nb::arg("method_id"), nb::arg("atom_types"),
               nb::arg("bond_types"), nb::arg("id"));
    module.def("_parameter_values", &parameter_set_values, nb::arg("catalog"), nb::arg("id"));
    module.def("_with_values", &parameter_set_with_values, nb::arg("catalog"), nb::arg("id"),
               nb::arg("values"), nb::arg("new_id"), nb::arg("name"), nb::arg("publication"),
               nb::arg("notes"));
}

} // namespace chargefw::python
