#pragma once

#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>

namespace chargefw::calculation {

// Which part of a request for an assessment or a calculation is invalid. More kinds may be added,
// so a switch needs a default branch.
enum class RequestErrorKind : std::uint8_t {
    // The method ID is not registered.
    unknown_method,
    // The execution mode, radius, or atom thresholds do not go together or are out of range.
    invalid_execution,
};

[[nodiscard]] auto to_string(RequestErrorKind value) -> std::string_view;

// A request that is invalid before any method or parameter set is considered, so that no choice
// of them could satisfy it. A valid request that a method or parameter set cannot satisfy is not
// an error: the assessment reports it as a rejection. kind() says which part of the request to
// change and reason(), which is also what(), why. option() is the ID of the method option the
// error concerns, when it concerns one.
class RequestError : public std::invalid_argument {
  public:
    RequestError(RequestErrorKind kind, const std::string& reason,
                 std::optional<std::string> option = {});

    [[nodiscard]] auto kind() const noexcept -> RequestErrorKind;
    [[nodiscard]] auto reason() const noexcept -> std::string_view;
    [[nodiscard]] auto option() const noexcept -> const std::optional<std::string>&;

  private:
    RequestErrorKind kind_;
    std::optional<std::string> option_;
};

} // namespace chargefw::calculation
