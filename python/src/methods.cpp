#include "bindings.h"

#include <chargefw/methods/method_prerequisites.h>
#include <chargefw/methods/method_registry.h>

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/variant.h>

#include <string>

namespace nb = nanobind;

namespace chargefw::python {
namespace {

auto method_descriptors() -> nb::list {
    auto result = nb::list{};
    for (const auto& method : methods::method_registry().methods()) {
        const auto& metadata = method->metadata();
        const auto requirements = method->requirements();
        auto descriptor = nb::dict{};
        descriptor["id"] = std::string{metadata.id};
        descriptor["name"] = std::string{metadata.name};
        descriptor["full_name"] = std::string{metadata.full_name};
        descriptor["publication"] = metadata.publication.has_value()
                                        ? nb::cast(std::string{*metadata.publication})
                                        : nb::none();
        descriptor["notes"] = std::string{metadata.notes};
        descriptor["recommended"] = metadata.recommended;
        descriptor["requires_coordinates"] = requirements.coordinates;
        descriptor["requires_bonds"] = requirements.bonds;
        descriptor["time_complexity"] =
            std::string{methods::to_string(requirements.resources.time)};
        descriptor["memory_complexity"] =
            std::string{methods::to_string(requirements.resources.memory)};
        descriptor["supports_cutoff"] = requirements.resources.supports_cutoff;
        descriptor["supports_cover"] = requirements.resources.supports_cover;
        descriptor["supports_fixed_ions"] = requirements.supports_fixed_point_sources;
        auto options = nb::list{};
        for (const auto& spec : method->option_schema()) {
            auto option = nb::dict{};
            option["id"] = std::string{spec.id};
            option["description"] = std::string{spec.description};
            option["type"] = std::string{methods::to_string(spec.type)};
            option["default"] = nb::cast(spec.default_value);
            auto choices = nb::list{};
            for (const auto& choice : spec.choices) {
                choices.append(nb::cast(choice));
            }
            option["choices"] = std::move(choices);
            option["minimum"] = spec.minimum.has_value() ? nb::cast(*spec.minimum) : nb::none();
            option["minimum_inclusive"] = spec.minimum_inclusive;
            option["maximum"] = spec.maximum.has_value() ? nb::cast(*spec.maximum) : nb::none();
            option["maximum_inclusive"] = spec.maximum_inclusive;
            options.append(std::move(option));
        }
        descriptor["options"] = std::move(options);
        result.append(std::move(descriptor));
    }
    return result;
}

} // namespace

void bind_methods(nb::module_& module) {
    module.def("_method_descriptors", &method_descriptors);
}

} // namespace chargefw::python
