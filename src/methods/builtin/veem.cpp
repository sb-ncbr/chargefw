#include "methods/builtin/veem.h"

#include "methods/builtin/prerequisite_helpers.h"

#include <chargefw/core/molecule.h>
#include <chargefw/core/periodic_table.h>
#include <chargefw/methods/method_failure.h>

#include <stdexcept>
#include <vector>

namespace chargefw::methods::builtin {

auto VEEMMethod::add_method_specific_prerequisite_issues(const MethodPrerequisiteInput& input,
                                                         PrerequisiteResult& result) const -> void {
    detail::add_element_prerequisite_issues(
        input, result, "the method requires valence electrons and electronegativity",
        [](const core::Element& element) -> bool {
            return element.valence_electron_count().has_value() && element.electronegativity != 0.0;
        });

    if (core::total_formal_charge(input.prepared_molecule.molecule()) != 0) {
        result.issues.push_back(PrerequisiteIssue{
            .kind = PrerequisiteIssueKind::unsupported_molecule,
            .message =
                "the method supports only neutral molecules because its valence-electron equation "
                "conserves zero total charge"});
    }
}

auto VEEMMethod::calculate(const CalculationInput& input) const -> charges::AtomicCharges {
    const auto& table = core::periodic_table();
    const auto& molecule = input.molecule();
    std::vector charges(molecule.atom_count(), 0.0);

    if (molecule.atom_count() == 0) {
        return charges::AtomicCharges{std::move(charges)};
    }

    auto numerator = 0.0;
    auto denominator = 0.0;

    for (const auto& atom : molecule.atoms()) {
        const auto& element = table.element(atom.atomic_number());
        const auto valence = static_cast<double>(*element.valence_electron_count());

        numerator += element.electronegativity * valence;
        denominator += valence;
    }

    if (denominator == 0.0) {
        throw MethodFailure{FailureKind::no_solution,
                            "the atoms have no valence electrons in total; use another method"};
    }

    const auto equalized_electronegativity = numerator / denominator;

    if (equalized_electronegativity == 0.0) {
        throw MethodFailure{FailureKind::no_solution,
                            "the equalized electronegativity is zero; use another method"};
    }

    for (std::size_t atom_index = 0; atom_index < molecule.atom_count(); ++atom_index) {
        const auto& atom = molecule.atom(atom_index);
        const auto& element = table.element(atom.atomic_number());
        const auto valence = static_cast<double>(*element.valence_electron_count());

        charges[atom_index] = valence * (equalized_electronegativity - element.electronegativity) /
                              equalized_electronegativity;
    }

    return charges::AtomicCharges{std::move(charges)};
}

auto VEEMMethod::estimate(const TargetSize& target, const MethodOptions& /*options*/) const
    -> ResourceEstimate {
    return {.memory_bytes = 0.0, .operations = target.atoms};
}

} // namespace chargefw::methods::builtin
