#pragma once

#include <nanobind/nanobind.h>

namespace chargefw::python {

void bind_core(nanobind::module_& module);
void bind_adapters(nanobind::module_& module);
void bind_methods(nanobind::module_& module);
void bind_parameters(nanobind::module_& module);
void bind_calculation(nanobind::module_& module);
void bind_parameterisation(nanobind::module_& module);

} // namespace chargefw::python
