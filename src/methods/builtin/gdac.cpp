#include "methods/builtin/gdac.h"

#include "methods/builtin/prerequisite_helpers.h"

#include <chargefw/parameters/models/parameter_view.h>

#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

namespace chargefw::methods::builtin {
namespace {

[[nodiscard]] auto vdw_radius(const core::Atom& atom) -> double {
    return core::periodic_table().element(atom.atomic_number()).van_der_waals_radius;
}

} // namespace

auto GDACMethod::add_method_specific_prerequisite_issues(const MethodPrerequisiteInput& input,
                                                         PrerequisiteResult& result) const -> void {
    detail::add_component_neutrality_prerequisite_issue(input, result, "GDAC");

    detail::add_element_prerequisite_issues(
        input, result, "GDAC requires a positive van der Waals radius",
        [](const core::Element& element) -> bool { return element.van_der_waals_radius > 0.0; });
}

auto GDACMethod::calculate(const CalculationInput& input) const -> charges::AtomicCharges {
    const auto iterations = input.method_options().get<int>("iters");

    const auto& molecule = input.molecule();
    const auto& geometry = input.geometry();
    const auto& parameters = input.parameters();

    const auto parameter_a = parameters.atom("A");
    const auto parameter_b = parameters.atom("B");

    const auto atom_count = molecule.atom_count();

    std::vector<double> charges(atom_count, 0.0);
    std::vector<double> electronegativities(atom_count, 0.0);

    for (int iteration = 0; iteration < iterations; ++iteration) {
        const auto alpha = iteration + 1;
        for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
            electronegativities[atom_index] =
                parameter_b[atom_index] * charges[atom_index] + parameter_a[atom_index];
        }

        for (const auto& bond : molecule.bonds()) {
            auto lower_atom_index = bond.first_atom_index();
            auto higher_atom_index = bond.second_atom_index();

            auto lower_electronegativity = electronegativities[lower_atom_index];
            auto higher_electronegativity = electronegativities[higher_atom_index];

            if (lower_electronegativity > higher_electronegativity) {
                std::swap(lower_atom_index, higher_atom_index);
                std::swap(lower_electronegativity, higher_electronegativity);
            }

            const auto& lower_atom = molecule.atom(lower_atom_index);
            const auto& higher_atom = molecule.atom(higher_atom_index);

            const auto denominator = parameter_a[lower_atom_index] + parameter_b[lower_atom_index];

            if (denominator == 0.0) {
                throw std::logic_error{"GDAC denominator is zero"};
            }

            const auto radius_sum = vdw_radius(lower_atom) + vdw_radius(higher_atom);

            const auto factor =
                1.0 - (geometry.distance(lower_atom_index, higher_atom_index) / radius_sum);

            const auto difference = std::pow(factor, alpha) *
                                    (higher_electronegativity - lower_electronegativity) /
                                    denominator;

            charges[lower_atom_index] += difference;
            charges[higher_atom_index] -= difference;
        }
    }

    return charges::AtomicCharges{std::move(charges)};
}

auto GDACMethod::estimate(const TargetSize& target, const MethodOptions& options) const
    -> ResourceEstimate {
    const auto iterations = options.get<int>("iters");
    // Electronegativity vector; each iteration visits every atom and bond.
    return {.memory_bytes = 8.0 * target.atoms,
            .operations = iterations * (target.atoms + target.bonds)};
}

} // namespace chargefw::methods::builtin
