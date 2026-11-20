#include <chargefw/io/parse_error.h>

#include "parse_context.h"

#include <istream>
#include <stdexcept>
#include <utility>

namespace chargefw::io {
namespace {

[[nodiscard]] auto describe(const std::string& reason, const ParseLocation& location)
    -> std::string {
    auto text = std::string{};
    if (location.line.has_value()) {
        text = (location.source_name.empty() ? "line " : location.source_name + ":") +
               std::to_string(*location.line) + ": ";
    } else if (!location.source_name.empty()) {
        text = location.source_name + ": ";
    }
    if (location.record_index.has_value()) {
        text +=
            "record " + std::to_string(*location.record_index + 1) +
            (location.record_id.empty() ? "" : " ('" + location.record_id.display_string() + "')") +
            ": ";
    }
    return text + reason;
}

} // namespace

auto to_string(const ParseErrorKind value) -> std::string_view {
    switch (value) {
    case ParseErrorKind::no_records:
        return "no_records";
    case ParseErrorKind::invalid_encoding:
        return "invalid_encoding";
    case ParseErrorKind::malformed:
        return "malformed";
    case ParseErrorKind::unsupported:
        return "unsupported";
    case ParseErrorKind::inconsistent_models:
        return "inconsistent_models";
    }
    throw std::logic_error{"unknown parse error kind"};
}

ParseError::ParseError(const ParseErrorKind kind, std::string reason, ParseLocation location)
    : std::runtime_error{describe(reason, location)}, kind_{kind}, reason_{std::move(reason)},
      location_{std::move(location)} {}

auto ParseError::kind() const noexcept -> ParseErrorKind {
    return kind_;
}

auto ParseError::reason() const noexcept -> const std::string& {
    return reason_;
}

auto ParseError::location() const noexcept -> const ParseLocation& {
    return location_;
}

namespace detail {

auto rethrow_located(const ParseLocation& location) -> void {
    try {
        throw;
    } catch (const ParseError& error) {
        auto known = error.location();
        if (known.source_name.empty()) {
            known.source_name = location.source_name;
        }
        if (!known.line.has_value()) {
            known.line = location.line;
        }
        if (!known.record_index.has_value()) {
            known.record_index = location.record_index;
        }
        if (known.record_id.empty()) {
            known.record_id = location.record_id;
        }
        throw ParseError{error.kind(), error.reason(), std::move(known)};
    } catch (const std::invalid_argument& error) {
        throw ParseError{ParseErrorKind::malformed, error.what(), location};
    }
}

auto require_readable(std::istream& input) -> void {
    if (input.bad() || (input.fail() && !input.eof())) {
        throw std::runtime_error{"failed to read the input"};
    }
}

} // namespace detail

} // namespace chargefw::io
