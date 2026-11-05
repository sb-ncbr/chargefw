#include "bindings.h"
#include "native_parameter_catalog.h"

#include <chargefw/core/molecule.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/core/periodic_table.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method.h>
#include <chargefw/methods/method_registry.h>
#include <chargefw/parameters/classification/parameter_template.h>
#include <chargefw/parameters/io/parameter_set_io.h>
#include <chargefw/parameters/models/parameter_key.h>
#include <chargefw/parameters/models/parameter_set_metadata.h>
#include <chargefw/parameters/models/parameter_values.h>

#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>

#include <cstddef>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace nb = nanobind;

namespace chargefw::python {
namespace {

// Parameter values passed from Python, one per slot in layout order.
using value_array = nb::ndarray<const double, nb::ndim<1>, nb::c_contig, nb::device::cpu>;

auto make_bundled_parameter_catalog() -> NativeParameterCatalog {
    return NativeParameterCatalog{parameters::load_default_parameter_sets()};
}

// A catalog with the one parameter set of a JSON file, for sets that are not bundled.
auto load_parameter_set_catalog(const std::string& path) -> NativeParameterCatalog {
    auto parameter_sets = std::vector<parameters::ParameterSet>{};
    parameter_sets.push_back(parameters::load_parameter_set_json_file(path));
    return NativeParameterCatalog{std::move(parameter_sets)};
}

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
    parameter_sets.push_back(parameters::make_parameter_template(
        prepared,
        {.common = parameter_names(requirements.common_parameters),
         .atom = parameter_names(requirements.atom_parameters),
         .bond = parameter_names(requirements.bond_parameters)},
        parameters::atom_classification_kind_from_string(atom_types),
        parameters::bond_classification_kind_from_string(bond_types), std::move(metadata)));
    return NativeParameterCatalog{std::move(parameter_sets)};
}

auto find_parameter_set(const NativeParameterCatalog& catalog, const std::string& id)
    -> const parameters::ParameterSet& {
    for (const auto& parameter_set : catalog.parameter_sets()) {
        if (parameter_set.id() == id) {
            return parameter_set;
        }
    }
    throw std::invalid_argument{"no parameter set '" + id + "' in the catalog"};
}

// The JSON text of the catalog's parameter set with the given ID.
auto parameter_set_json(const NativeParameterCatalog& catalog, const std::string& id)
    -> std::string {
    auto output = std::ostringstream{};
    parameters::save_parameter_set_json(find_parameter_set(catalog, id), output);
    return output.str();
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
    parameter_sets.push_back(parameters::with_parameter_values(
        find_parameter_set(catalog, id), std::span{values.data(), values.size()},
        std::move(metadata)));
    return NativeParameterCatalog{std::move(parameter_sets)};
}

auto atom_key_payload(const parameters::AtomParameterKey& key) -> nb::dict {
    auto payload = nb::dict{};
    payload["element"] = key.atomic_number == 0
                             ? std::string{"*"}
                             : std::string{core::element_symbol(key.atomic_number)};
    payload["classifier"] = std::string{parameters::to_string(key.classification)};
    payload["type"] = key.type;
    return payload;
}

auto parameter_set_descriptors(const NativeParameterCatalog& catalog) -> nb::list {
    auto result = nb::list{};
    for (const auto& parameter_set : catalog.parameter_sets()) {
        const auto& metadata = parameter_set.metadata();
        auto descriptor = nb::dict{};
        descriptor["id"] = metadata.id;
        descriptor["method_id"] = metadata.method_id;
        descriptor["name"] = metadata.name;
        descriptor["publication"] = metadata.publication;
        descriptor["notes"] = metadata.notes;
        auto atom_keys = nb::list{};
        for (const auto& entry : parameter_set.atom().entries()) {
            atom_keys.append(atom_key_payload(entry.key));
        }
        descriptor["atom_keys"] = std::move(atom_keys);
        auto bond_keys = nb::list{};
        for (const auto& entry : parameter_set.bond().entries()) {
            auto payload = nb::dict{};
            payload["atoms"] = nb::make_tuple(atom_key_payload(entry.key.first_atom),
                                              atom_key_payload(entry.key.second_atom));
            payload["classifier"] =
                std::string{parameters::to_string(entry.key.bond.classification)};
            payload["type"] = entry.key.bond.type;
            bond_keys.append(std::move(payload));
        }
        descriptor["bond_keys"] = std::move(bond_keys);
        result.append(std::move(descriptor));
    }
    return result;
}

} // namespace

void bind_parameters(nb::module_& module) {
    nb::class_<NativeParameterCatalog>(module, "_NativeParameterCatalog")
        .def("_descriptors", &parameter_set_descriptors)
        .def("_json", &parameter_set_json, nb::arg("id"))
        .def("_with_values", &parameter_set_with_values, nb::arg("id"), nb::arg("values"),
             nb::arg("new_id"), nb::arg("name"), nb::arg("publication"), nb::arg("notes"));
    module.def("_load_bundled_parameter_catalog", &make_bundled_parameter_catalog,
               nb::call_guard<nb::gil_scoped_release>());
    module.def("_load_parameter_set_catalog", &load_parameter_set_catalog, nb::arg("path"),
               nb::call_guard<nb::gil_scoped_release>());
    module.def("_make_parameter_template_catalog", &make_parameter_template_catalog,
               nb::arg("molecules"), nb::arg("method_id"), nb::arg("atom_types"),
               nb::arg("bond_types"), nb::arg("id"));
}

} // namespace chargefw::python
