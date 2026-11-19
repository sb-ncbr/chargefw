#pragma once

#include <chargefw/core/issue_location.h>
#include <chargefw/parameters/classification/parameter_classification.h>

#include <cstddef>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace chargefw::parameters {

enum class ClassificationIssueKind : std::uint8_t {
    missing_atom_parameter,
    missing_bond_parameter,
};

// An atom or bond of the classified molecule that no entry matches. The message names it but not
// the parameter set; `location` has its atom_index or bond_index.
struct ClassificationIssue {
    ClassificationIssueKind kind;
    std::string message;
    core::IssueLocation location = {};
};

// Strict matching uses only entries for an atom's or bond's own type; permissive matching then
// tries the next lower bond order.
enum class ParameterMatching : std::uint8_t {
    strict,
    permissive,
};

[[nodiscard]] auto parameter_matching_from_string(std::string_view value) -> ParameterMatching;
[[nodiscard]] auto to_string(ParameterMatching value) -> std::string_view;

struct ClassificationOptions {
    ParameterMatching parameter_matching = ParameterMatching::strict;
};

class ClassificationResult {
  public:
    explicit ClassificationResult(ParameterClassification classification);
    explicit ClassificationResult(std::vector<ClassificationIssue> issues);

    [[nodiscard]] explicit operator bool() const noexcept;

    [[nodiscard]] auto classification() const -> const ParameterClassification&;
    [[nodiscard]] auto issues() const noexcept -> const std::vector<ClassificationIssue>&;

  private:
    ParameterClassification classification_;
    std::vector<ClassificationIssue> issues_;
};

} // namespace chargefw::parameters
