#include "methods/builtin/eqeqc.h"
#include "methods/atom_parameter_values.h"

#include "methods/builtin/eqeq.h"

#include <chargefw/core/molecule.h>
#include <chargefw/core/periodic_table.h>
#include <chargefw/parameters/models/parameter_view.h>

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace chargefw::methods::builtin {
namespace {

[[nodiscard]] auto correction_term(const double alpha, const double distance,
                                   const double first_covalent_radius,
                                   const double second_covalent_radius, const double first_dz,
                                   const double second_dz) -> double {
    const auto tkk = first_dz - second_dz;
    const auto bkk = std::exp(-alpha * (distance - first_covalent_radius - second_covalent_radius));

    return tkk * bkk;
}

} // namespace

auto EQeqCMethod::calculate(const CalculationInput& input) const -> charges::AtomicCharges {
    const auto& molecule = input.molecule();
    const auto& geometry = input.geometry();
    const auto& parameters = input.parameters();

    const auto atom_count = molecule.atom_count();

    if (atom_count == 0) {
        return charges::AtomicCharges{std::vector<double>{}};
    }

    const auto alpha = parameters.common("alpha");

    if (alpha <= 0.0) {
        throw std::logic_error{"EQeq+C common parameter 'alpha' must be positive"};
    }

    // The pair loop below reads Dz O(n^2) times; resolve it once per atom.
    const auto dz = methods::detail::atom_parameter_values(parameters, "Dz", atom_count);
    auto values = eqeq_core::calculate(input, "EQeq+C");

    const auto& table = core::periodic_table();

    for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
        const auto& atom = molecule.atom(atom_index);
        const auto& element = table.element(atom.atomic_number());

        auto correction = 0.0;

        for (std::size_t other_atom_index = 0; other_atom_index < atom_count; ++other_atom_index) {
            if (atom_index == other_atom_index) {
                continue;
            }

            const auto& other_atom = molecule.atom(other_atom_index);
            const auto& other_element = table.element(other_atom.atomic_number());
            const auto distance = geometry.distance(atom_index, other_atom_index);

            correction += correction_term(alpha, distance, element.covalent_radius,
                                          other_element.covalent_radius, dz[atom_index],
                                          dz[other_atom_index]);
        }

        values[atom_index] += correction;
    }

    return charges::AtomicCharges{std::move(values)};
}

auto EQeqCMethod::estimate(const TargetSize& target, const MethodOptions& /*options*/) const
    -> ResourceEstimate {
    const auto n = target.atoms;
    const auto size = n + 1.0;
    // Dense (n + 1) x (n + 1) system with a charge constraint; the correction adds an all-pairs
    // pass.
    return {.memory_bytes = 8.0 * size * size,
            .operations = 2.0 / 3.0 * size * size * size + n * n + 0.5 * n * n};
}

} // namespace chargefw::methods::builtin
