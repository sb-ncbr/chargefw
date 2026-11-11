#include <chargefw/features/prepared_molecule.h>
#include <chargefw/methods/method.h>
#include <chargefw/methods/parameter_prerequisites.h>
#include <chargefw/parameters/classification/parameter_classifier.h>

#include "core/diagnostic_description.h"
#include "methods/parameter_issue_messages.h"

#include <cstddef>
#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace chargefw::methods {
namespace {

struct ParameterPrerequisiteResult {
    std::vector<PrerequisiteIssue> issues;
    std::optional<parameters::ParameterClassification> classification;

    [[nodiscard]] explicit operator bool() const noexcept {
        return issues.empty();
    }
};

auto add_missing_parameter_issue(ParameterPrerequisiteResult& result, std::string message) -> void {
    result.issues.push_back(PrerequisiteIssue{.kind = PrerequisiteIssueKind::missing_parameters,
                                              .message = std::move(message)});
}

auto check_common_parameters(const MethodRequirements& requirements,
                             const parameters::ParameterSet& parameter_set,
                             ParameterPrerequisiteResult& result) -> void {
    for (const auto name : requirements.common_parameters) {
        if (!parameter_set.common().contains(name)) {
            add_missing_parameter_issue(
                result, "parameter set '" + std::string{parameter_set.id()} +
                            "' is missing common parameter '" + std::string{name} + "'");
        }
    }
}

auto check_atom_parameter_names(const MethodRequirements& requirements,
                                const parameters::ParameterSet& parameter_set,
                                ParameterPrerequisiteResult& result) -> void {
    if (!requirements.requires_atom_parameters()) {
        return;
    }

    if (parameter_set.atom().empty()) {
        add_missing_parameter_issue(result, "parameter set '" + std::string{parameter_set.id()} +
                                                "' has no atom parameter entries");
        return;
    }

    for (std::size_t entry_index = 0; entry_index < parameter_set.atom().size(); ++entry_index) {
        for (const auto name : requirements.atom_parameters) {
            if (!parameter_set.atom().contains(entry_index, name)) {
                add_missing_parameter_issue(
                    result, "atom parameter entry " + std::to_string(entry_index) +
                                " in parameter set '" + std::string{parameter_set.id()} +
                                "' is missing parameter '" + std::string{name} + "'");
            }
        }
    }
}

auto check_bond_parameter_names(const MethodRequirements& requirements,
                                const parameters::ParameterSet& parameter_set,
                                ParameterPrerequisiteResult& result) -> void {
    if (!requirements.requires_bond_parameters()) {
        return;
    }

    if (parameter_set.bond().empty()) {
        add_missing_parameter_issue(result, "parameter set '" + std::string{parameter_set.id()} +
                                                "' has no bond parameter entries");
        return;
    }

    for (std::size_t entry_index = 0; entry_index < parameter_set.bond().size(); ++entry_index) {
        for (const auto name : requirements.bond_parameters) {
            if (!parameter_set.bond().contains(entry_index, name)) {
                add_missing_parameter_issue(
                    result, "bond parameter entry " + std::to_string(entry_index) +
                                " in parameter set '" + std::string{parameter_set.id()} +
                                "' is missing parameter '" + std::string{name} + "'");
            }
        }
    }
}

// One issue for the atoms without parameters and one for the bonds, not one per atom or bond.
auto add_classification_issues(ParameterPrerequisiteResult& result,
                               const std::vector<parameters::ClassificationIssue>& issues,
                               const core::Molecule& molecule,
                               const std::string_view parameter_set_id) -> void {
    auto atoms = std::vector<std::size_t>{};
    auto bonds = std::vector<std::size_t>{};
    for (const auto& issue : issues) {
        switch (issue.kind) {
        case parameters::ClassificationIssueKind::MISSING_ATOM_PARAMETER:
            atoms.push_back(issue.object_index);
            break;
        case parameters::ClassificationIssueKind::MISSING_BOND_PARAMETER:
            bonds.push_back(issue.object_index);
            break;
        }
    }
    if (!atoms.empty()) {
        result.issues.push_back(PrerequisiteIssue{
            .kind = PrerequisiteIssueKind::parameter_classification_failed,
            .message = detail::missing_atom_parameters_message(molecule, parameter_set_id, atoms),
            .location = {.atom_index = atoms.front(), .atom_indices = std::move(atoms)}});
    }
    if (!bonds.empty()) {
        result.issues.push_back(PrerequisiteIssue{
            .kind = PrerequisiteIssueKind::parameter_classification_failed,
            .message = detail::missing_bond_parameters_message(molecule, parameter_set_id, bonds),
            .location = {.bond_index = bonds.front(), .bond_indices = std::move(bonds)}});
    }
}

auto check_molecule_parameters(const Method& method,
                               const features::PreparedMolecule& prepared_molecule,
                               const parameters::ParameterSet& parameter_set,
                               const parameters::ClassificationOptions& classification_options)
    -> ParameterPrerequisiteResult {
    ParameterPrerequisiteResult result;

    const auto requirements = method.requirements();

    if (!requirements.requires_parameters()) {
        return result;
    }

    if (!parameter_set.method_id().empty() && parameter_set.method_id() != method.id()) {
        add_missing_parameter_issue(result, "parameter set '" + std::string{parameter_set.id()} +
                                                "' belongs to method '" +
                                                std::string{parameter_set.method_id()} +
                                                "', not '" + std::string{method.id()} + "'");
        return result;
    }

    check_common_parameters(requirements, parameter_set, result);
    check_atom_parameter_names(requirements, parameter_set, result);
    check_bond_parameter_names(requirements, parameter_set, result);

    if (!result.issues.empty()) {
        return result;
    }

    if (requirements.requires_atom_parameters() || requirements.requires_bond_parameters()) {
        const auto classification = parameters::try_classify_parameters(
            prepared_molecule.topology(), parameter_set, classification_options);

        if (!classification) {
            add_classification_issues(result, classification.issues(), prepared_molecule.molecule(),
                                      parameter_set.id());
            return result;
        }

        result.classification = classification.classification();
    }

    return result;
}

} // namespace

auto check_parameter_prerequisites(const Method& method,
                                   const features::PreparedMoleculeCollection& molecules,
                                   const parameters::ParameterSet& parameter_set,
                                   const parameters::ClassificationOptions& classification_options)
    -> CollectionParameterPrerequisiteResult {
    auto result = CollectionParameterPrerequisiteResult{};
    result.classifications.reserve(molecules.size());

    for (std::size_t molecule_index = 0; molecule_index < molecules.size(); ++molecule_index) {
        const auto molecule_result = check_molecule_parameters(
            method, molecules[molecule_index], parameter_set, classification_options);

        if (!molecule_result) {
            for (const auto& issue : molecule_result.issues) {
                result.issues.push_back(PrerequisiteIssue{
                    .kind = issue.kind,
                    .message = core::detail::molecule_description(
                                   molecules[molecule_index].molecule(), molecule_index) +
                               ": " + issue.message,
                    .molecule_index = molecule_index,
                    .location = issue.location});
            }

            continue;
        }

        if (molecule_result.classification.has_value()) {
            result.classifications.push_back(*molecule_result.classification);
        }
    }

    if (!result.issues.empty()) {
        result.classifications.clear();
    }

    return result;
}

} // namespace chargefw::methods
