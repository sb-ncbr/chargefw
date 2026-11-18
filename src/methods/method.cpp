#include <chargefw/methods/method.h>

#include "methods/geometry_issue_messages.h"

#include <chargefw/features/prepared_molecule.h>

#include <exception>
#include <stdexcept>
#include <string>
#include <utility>

namespace chargefw::methods {
namespace {

auto add_invalid_options_issue(PrerequisiteResult& result, const std::exception& error) -> void {
    result.issues.push_back(
        PrerequisiteIssue{.kind = PrerequisiteIssueKind::invalid_options, .message = error.what()});
}

auto add_coordinate_prerequisite_issues(const features::PreparedMolecule& prepared,
                                        PrerequisiteResult& result) -> void {
    const auto& molecule = prepared.molecule();
    const auto defects = prepared.geometry_defects();
    for (std::size_t conformer_index = 0; conformer_index < defects.size(); ++conformer_index) {
        const auto& nonfinite_atom = defects[conformer_index].nonfinite_atom;
        const auto& coincident = defects[conformer_index].coincident_atoms;

        if (nonfinite_atom.has_value()) {
            auto location = core::IssueLocation{.conformer_index = conformer_index,
                                                .atom_index = nonfinite_atom,
                                                .atom_indices = {*nonfinite_atom}};
            auto message = detail::invalid_geometry_message(molecule, location);
            result.issues.push_back(
                PrerequisiteIssue{.kind = PrerequisiteIssueKind::invalid_geometry,
                                  .message = std::move(message),
                                  .location = std::move(location)});
            continue;
        }

        if (coincident.has_value()) {
            auto location =
                core::IssueLocation{.conformer_index = conformer_index,
                                    .atom_index = coincident->second,
                                    .other_atom_index = coincident->first,
                                    .atom_indices = {coincident->first, coincident->second}};
            auto message = detail::invalid_geometry_message(molecule, location);
            result.issues.push_back(
                PrerequisiteIssue{.kind = PrerequisiteIssueKind::invalid_geometry,
                                  .message = std::move(message),
                                  .location = std::move(location)});
        }
    }
}

} // namespace

auto to_string(const ComplexityTerm value) -> std::string_view {
    switch (value) {
    case ComplexityTerm::constant:
        return "O(1)";
    case ComplexityTerm::atoms:
        return "O(n)";
    case ComplexityTerm::bonds:
        return "O(m)";
    case ComplexityTerm::atoms_plus_bonds:
        return "O(n + m)";
    case ComplexityTerm::atoms_squared:
        return "O(n^2)";
    case ComplexityTerm::atoms_cubed:
        return "O(n^3)";
    case ComplexityTerm::bonds_squared:
        return "O(m^2)";
    case ComplexityTerm::bonds_cubed:
        return "O(m^3)";
    case ComplexityTerm::atoms_plus_bonds_squared:
        return "O((n + m)^2)";
    case ComplexityTerm::atoms_plus_bonds_cubed:
        return "O((n + m)^3)";
    case ComplexityTerm::atoms_squared_plus_bonds_squared:
        return "O(n^2 + m^2)";
    }
    throw std::invalid_argument{"unknown complexity term"};
}

auto Method::requires_parameters() const -> bool {
    return requirements().requires_parameters();
}

auto Method::check_method_prerequisites(const MethodPrerequisiteInput& input) const
    -> PrerequisiteResult {
    PrerequisiteResult result;

    try {
        validate_method_options(option_schema(), input.method_options);
    } catch (const std::exception& error) {
        add_invalid_options_issue(result, error);
    }

    const auto method_requirements = requirements();
    const auto& molecule = input.prepared_molecule.molecule();

    if (method_requirements.coordinates && !molecule.has_coordinates()) {
        result.issues.push_back(PrerequisiteIssue{
            .kind = PrerequisiteIssueKind::missing_feature,
            .message = "the method requires coordinates, but the molecule has none"});
    }
    if (method_requirements.coordinates && molecule.has_coordinates()) {
        add_coordinate_prerequisite_issues(input.prepared_molecule, result);
    }
    if (method_requirements.bonds && molecule.atom_count() > 1 && molecule.bond_count() == 0) {
        result.issues.push_back(
            PrerequisiteIssue{.kind = PrerequisiteIssueKind::missing_feature,
                              .message = "the method requires bonds, but the molecule has none"});
    }

    add_method_specific_prerequisite_issues(input, result);

    return result;
}

auto Method::add_method_specific_prerequisite_issues(const MethodPrerequisiteInput& /* unused */,
                                                     PrerequisiteResult& /* unused */) const
    -> void {}

} // namespace chargefw::methods
