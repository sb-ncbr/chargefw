#include <chargefw/methods/method.h>

#include "core/diagnostic_description.h"

#include <chargefw/features/prepared_molecule.h>

#include <exception>
#include <stdexcept>
#include <string>

namespace chargefw::methods {
namespace {

auto add_invalid_options_issue(PrerequisiteResult& result, const std::exception& error) -> void {
    result.issues.push_back(
        PrerequisiteIssue{.kind = PrerequisiteIssueKind::invalid_options, .message = error.what()});
}

auto add_coordinate_prerequisite_issues(const Method& method,
                                        const features::PreparedMolecule& prepared,
                                        PrerequisiteResult& result) -> void {
    const auto& molecule = prepared.molecule();
    const auto defects = prepared.geometry_defects();
    for (std::size_t conformer_index = 0; conformer_index < defects.size(); ++conformer_index) {
        const auto& nonfinite_atom = defects[conformer_index].nonfinite_atom;
        const auto& coincident = defects[conformer_index].coincident_atoms;

        if (nonfinite_atom.has_value()) {
            result.issues.push_back(PrerequisiteIssue{
                .kind = PrerequisiteIssueKind::invalid_geometry,
                .message = "method '" + std::string{method.id()} + "', conformer " +
                           std::to_string(conformer_index + 1) + ": " +
                           core::detail::atom_description(molecule, *nonfinite_atom) +
                           " has non-finite coordinates",
                .location = {.conformer_index = conformer_index,
                             .atom_index = nonfinite_atom,
                             .atom_indices = {*nonfinite_atom}}});
            continue;
        }

        if (coincident.has_value()) {
            result.issues.push_back(PrerequisiteIssue{
                .kind = PrerequisiteIssueKind::invalid_geometry,
                .message = "method '" + std::string{method.id()} + "', conformer " +
                           std::to_string(conformer_index + 1) + ": " +
                           core::detail::atom_description(molecule, coincident->first) + " and " +
                           core::detail::atom_description(molecule, coincident->second) +
                           " have identical coordinates",
                .location = {.conformer_index = conformer_index,
                             .atom_index = coincident->second,
                             .other_atom_index = coincident->first,
                             .atom_indices = {coincident->first, coincident->second}}});
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
        result.issues.push_back(PrerequisiteIssue{.kind = PrerequisiteIssueKind::missing_feature,
                                                  .message = "method '" + std::string{id()} +
                                                             "' requires coordinates"});
    }
    if (method_requirements.coordinates && molecule.has_coordinates()) {
        add_coordinate_prerequisite_issues(*this, input.prepared_molecule, result);
    }
    if (method_requirements.bonds && molecule.atom_count() > 1 && molecule.bond_count() == 0) {
        result.issues.push_back(
            PrerequisiteIssue{.kind = PrerequisiteIssueKind::missing_feature,
                              .message = "method '" + std::string{id()} +
                                         "' requires bonds, but the molecule has none"});
    }

    add_method_specific_prerequisite_issues(input, result);

    return result;
}

auto Method::add_method_specific_prerequisite_issues(const MethodPrerequisiteInput& /* unused */,
                                                     PrerequisiteResult& /* unused */) const
    -> void {}

} // namespace chargefw::methods
