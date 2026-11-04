#pragma once

#include "core/diagnostic_description.h"

#include "features/topology_helpers.h"

#include <chargefw/core/molecule.h>
#include <chargefw/core/periodic_table.h>
#include <chargefw/methods/method_prerequisites.h>

#include <concepts>
#include <cstddef>
#include <numeric>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace chargefw::methods::builtin::detail {

// Rejects molecules with a formally charged connected component. A non-empty reason is appended
// to the message as an explanation.
inline auto add_component_neutrality_prerequisite_issue(const MethodPrerequisiteInput& input,
                                                        PrerequisiteResult& result,
                                                        const std::string_view method_name,
                                                        const std::string_view reason = {})
    -> void {
    const auto& molecule = input.prepared_molecule.molecule();
    const auto components =
        features::connected_components(input.prepared_molecule.topology().adjacency());

    for (const auto& component : components) {
        const auto formal_charge =
            std::accumulate(component.begin(), component.end(), 0,
                            [&molecule](const int sum, const std::size_t atom_index) {
                                return sum + molecule.atom(atom_index).formal_charge();
                            });
        if (formal_charge != 0) {
            result.add(PrerequisiteIssue{
                .kind = PrerequisiteIssueKind::unsupported_molecule,
                .message = std::string{method_name} +
                           " supports only molecules in which every bonded part has zero total "
                           "formal charge" +
                           (reason.empty() ? std::string{} : " because " + std::string{reason})});
            return;
        }
    }
}

template <typename ElementSupported>
    requires std::predicate<ElementSupported, const core::Element&>

auto add_element_prerequisite_issues(const MethodPrerequisiteInput& input,
                                     PrerequisiteResult& result, const std::string_view requirement,
                                     const ElementSupported& element_supported) -> void {
    const auto& table = core::periodic_table();
    const auto& molecule = input.prepared_molecule.molecule();

    auto atoms = std::vector<std::size_t>{};
    for (std::size_t atom_index = 0; atom_index < molecule.atom_count(); ++atom_index) {
        const auto atomic_number = molecule.atom(atom_index).atomic_number();
        if (!table.contains(atomic_number) || !element_supported(table.element(atomic_number))) {
            atoms.push_back(atom_index);
        }
    }
    if (atoms.empty()) {
        return;
    }

    // One issue for all atoms whose element the method does not support.
    const auto example = core::detail::atom_description(molecule, atoms.front());
    result.add(PrerequisiteIssue{
        .kind = PrerequisiteIssueKind::unsupported_molecule,
        .message = atoms.size() == 1
                       ? std::string{requirement} + ": " + example
                       : std::string{requirement} + ": " + std::to_string(atoms.size()) + " of " +
                             std::to_string(molecule.atom_count()) + " atoms do not meet it (" +
                             core::detail::element_counts_description(molecule, atoms) +
                             "), such as " + example,
        .location = {.atom_index = atoms.front(), .atom_indices = std::move(atoms)}});
}

} // namespace chargefw::methods::builtin::detail
