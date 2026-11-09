#pragma once

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>

#include <memory>
#include <span>
#include <vector>

namespace chargefw::python {

// Values passed from Python as a contiguous one-dimensional array, such as parameter values in
// layout order.
using value_array =
    nanobind::ndarray<const double, nanobind::ndim<1>, nanobind::c_contig, nanobind::device::cpu>;

// A NumPy array that owns a copy of `values`.
[[nodiscard]] inline auto numpy_values(const std::span<const double> values)
    -> nanobind::ndarray<nanobind::numpy, const double, nanobind::shape<-1>> {
    auto owned_values = std::make_unique<std::vector<double>>(values.begin(), values.end());
    auto owner = nanobind::capsule{owned_values.get(), [](void* pointer) noexcept {
                                       delete static_cast<std::vector<double>*>(pointer);
                                   }};
    auto* released_values = owned_values.release();
    return {released_values->data(), {released_values->size()}, owner};
}

} // namespace chargefw::python
