#include "methods/builtin/peoe.h"

#include "methods/builtin/initial_charges.h"
#include "methods/builtin/prerequisite_helpers.h"

#include <chargefw/core/molecule.h>
#include <chargefw/parameters/models/parameter_view.h>

#include <cmath>
#include <utility>
#include <vector>

namespace chargefw::methods::builtin {

auto PEOEMethod::add_method_specific_prerequisite_issues(const MethodPrerequisiteInput& input,
                                                         PrerequisiteResult& result) const -> void {
    if (!detail::uses_formal_initial_charges(input.method_options)) {
        detail::add_component_neutrality_prerequisite_issue(input, result, "PEOE");
    }
}

auto PEOEMethod::calculate(const CalculationInput& input) const -> charges::AtomicCharges {
    const auto iterations = input.method_options().get<int>("iters");

    const auto& molecule = input.molecule();
    const auto& parameter_view = input.parameters();

    const auto parameter_a = parameter_view.atom("A");
    const auto parameter_b = parameter_view.atom("B");
    const auto parameter_c = parameter_view.atom("C");
    const auto damping_hydrogen = parameter_view.common("dampH");

    const auto atom_count = molecule.atom_count();

    auto charges = input.method_options().get<std::string>("initial_charges") == "formal"
                       ? detail::formal_initial_charges(molecule)
                       : std::vector<double>(atom_count, 0.0);
    std::vector electronegativities(atom_count, 0.0);

    for (int iteration = 0; iteration < iterations; ++iteration) {
        const auto alpha = iteration + 1;
        for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
            const auto charge = charges[atom_index];

            electronegativities[atom_index] = parameter_c[atom_index] * charge * charge +
                                              parameter_b[atom_index] * charge +
                                              parameter_a[atom_index];
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

            const auto denominator = lower_atom.atomic_number() == 1
                                         ? damping_hydrogen
                                         : parameter_a[lower_atom_index] +
                                               parameter_b[lower_atom_index] +
                                               parameter_c[lower_atom_index];

            const auto difference = std::pow(0.5, alpha) *
                                    (higher_electronegativity - lower_electronegativity) /
                                    denominator;

            charges[lower_atom_index] += difference;
            charges[higher_atom_index] -= difference;
        }
    }

    return charges::AtomicCharges{std::move(charges)};
}

auto PEOEMethod::estimate(const TargetSize& target, const MethodOptions& options) const
    -> ResourceEstimate {
    const auto iterations = options.get<int>("iters");
    // Electronegativity vector; each iteration visits every atom and bond.
    return {.memory_bytes = 8.0 * target.atoms,
            .operations = iterations * (target.atoms + target.bonds)};
}

} // namespace chargefw::methods::builtin
