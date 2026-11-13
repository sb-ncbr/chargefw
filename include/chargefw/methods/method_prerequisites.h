#pragma once

#include <chargefw/core/issue_location.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/methods/method_options.h>

#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace chargefw::methods {

class Method;

enum class PrerequisiteIssueKind : std::uint8_t {
    invalid_options,
    missing_feature,
    invalid_geometry,
    unsupported_molecule,
    missing_parameters,
    unknown_parameters,
    parameter_classification_failed,
    unsupported_fixed_ions,
};

[[nodiscard]] auto to_string(PrerequisiteIssueKind value) -> std::string_view;

struct PrerequisiteIssue {
    PrerequisiteIssueKind kind;
    std::string message;
    // Structured indices are zero-based. Human-readable messages use one-based numbering.
    std::optional<std::size_t> molecule_index = std::nullopt;
    core::IssueLocation location = {};
};

struct PrerequisiteResult {
    std::vector<PrerequisiteIssue> issues;

    [[nodiscard]] explicit operator bool() const noexcept {
        return issues.empty();
    }
};

struct MethodPrerequisiteInput {
    // NOLINTBEGIN(cppcoreguidelines-avoid-const-or-ref-data-members)
    const features::PreparedMolecule& prepared_molecule;
    const MethodOptions& method_options;
    // NOLINTEND(cppcoreguidelines-avoid-const-or-ref-data-members)
};

[[nodiscard]] auto check_method_prerequisites(const Method& method,
                                              const features::PreparedMoleculeCollection& molecules,
                                              const MethodOptions& method_options)
    -> PrerequisiteResult;

} // namespace chargefw::methods
