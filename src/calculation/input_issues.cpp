#include <chargefw/calculation/input_issues.h>

#include "input_checks.h"

#include <algorithm>
#include <stdexcept>

namespace chargefw::calculation {

auto to_string(const InputIssueKind value) -> std::string_view {
    switch (value) {
    case InputIssueKind::hydrogens_missing:
        return "hydrogens_missing";
    }
    throw std::invalid_argument{"unknown input issue kind"};
}

namespace detail {

auto find_input_issues(const core::MoleculeCollection& molecules) -> std::vector<InputIssue> {
    constexpr auto hydrogen = 1;
    constexpr auto carbon = 6;
    auto result = std::vector<InputIssue>{};
    for (std::size_t index = 0; index < molecules.size(); ++index) {
        const auto atoms = molecules[index].atoms();
        const auto has = [atoms](const int atomic_number) {
            return std::ranges::any_of(atoms, [atomic_number](const auto& atom) {
                return atom.atomic_number() == atomic_number;
            });
        };
        // Hydrogens are missing from a crystal structure before protonation or from an RDKit
        // molecule without AddHs(); a molecule such as CO2 or CCl4 has none, so this only warns.
        if (has(carbon) && !has(hydrogen)) {
            result.push_back(InputIssue{
                .kind = InputIssueKind::hydrogens_missing,
                .message = "the molecule contains carbon but no hydrogen atoms; if its hydrogens "
                           "were not added, add them before calculating",
                .molecule_index = index});
        }
    }
    return result;
}

} // namespace detail
} // namespace chargefw::calculation
