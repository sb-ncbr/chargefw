#pragma once

#include <chargefw/parameters/classification/parameter_classification.h>

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace chargefw::parameters {

enum class ClassificationIssueKind : std::uint8_t {
    missing_atom_parameter,
    missing_bond_parameter,
};

struct ClassificationIssue {
    ClassificationIssueKind kind;
    std::size_t object_index = 0;
    std::string message;
};

struct ClassificationOptions {
    bool permissive_types = false;
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
