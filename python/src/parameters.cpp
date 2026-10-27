#include "bindings.h"
#include "native_parameter_catalog.h"

#include <chargefw/core/periodic_table.h>
#include <chargefw/parameters/io/parameter_set_io.h>
#include <chargefw/parameters/models/parameter_key.h>

#include <nanobind/stl/string.h>

#include <string>

namespace nb = nanobind;

namespace chargefw::python {
namespace {

auto make_bundled_parameter_catalog() -> NativeParameterCatalog {
    return NativeParameterCatalog{parameters::load_default_parameter_sets()};
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
        descriptor["priority"] = metadata.priority;
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
        .def("_descriptors", &parameter_set_descriptors);
    module.def("_load_bundled_parameter_catalog", &make_bundled_parameter_catalog,
               nb::call_guard<nb::gil_scoped_release>());
}

} // namespace chargefw::python
