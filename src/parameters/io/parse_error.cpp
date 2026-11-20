#include <chargefw/parameters/io/parse_error.h>

#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

namespace chargefw::parameters {

auto to_string(const ParseErrorKind value) -> std::string_view {
    switch (value) {
    case ParseErrorKind::malformed:
        return "malformed";
    case ParseErrorKind::unsupported:
        return "unsupported";
    }
    throw std::logic_error{"unknown parameter parse error kind"};
}

ParseError::ParseError(const ParseErrorKind kind, std::string reason, std::string source_name)
    : std::runtime_error{source_name.empty() ? reason : source_name + ": " + reason}, kind_{kind},
      reason_{std::move(reason)}, source_name_{std::move(source_name)} {}

auto ParseError::kind() const noexcept -> ParseErrorKind {
    return kind_;
}

auto ParseError::reason() const noexcept -> const std::string& {
    return reason_;
}

auto ParseError::source_name() const noexcept -> const std::string& {
    return source_name_;
}

} // namespace chargefw::parameters
