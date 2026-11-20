#include "common_input.h"

#include "io/text_encoding.h"

#include <chargefw/io/parse_error.h>

#include <gemmi/atof.hpp>

#include <charconv>
#include <stdexcept>
#include <string>
#include <utility>

namespace chargefw::io::native::common_input {

auto rethrow_with_context(const MoleculeRecordIdentity& identity,
                          const std::optional<std::size_t> line) -> void {
    const auto location = ParseLocation{.source_name = identity.source_name,
                                        .line = line,
                                        .record_index = identity.record_index,
                                        .record_id = identity.record_id};
    try {
        throw;
    } catch (const ParseError& error) {
        // The throw site knows best where the offending data is, such as the line of a bond found
        // wrong at the end of its record; only what it left out is added.
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

auto at_end(std::istream& input) -> bool {
    require_readable(input);
    if (input.eof()) {
        return true;
    }
    const auto end = input.peek() == std::char_traits<char>::eof();
    require_readable(input);
    return end;
}

auto try_read_line(std::istream& input, std::size_t& line) -> std::optional<std::string> {
    std::string result;

    if (!std::getline(input, result)) {
        require_readable(input);
        return std::nullopt;
    }

    if (result.ends_with('\r')) {
        result.pop_back();
    }

    ++line;
    if (io::detail::first_non_ascii_line(result).has_value()) {
        throw ParseError{ParseErrorKind::invalid_encoding,
                         "non-ASCII character; MOL, SDF, and MOL2 files must be ASCII text",
                         {}};
    }
    return result;
}

auto read_line(std::istream& input, std::size_t& line, const std::string_view record_name)
    -> std::string {
    auto result = try_read_line(input, line);
    if (!result.has_value()) {
        ++line;
        throw ParseError{ParseErrorKind::malformed,
                         "unexpected end of " + std::string{record_name} + " record",
                         {}};
    }
    return std::move(*result);
}

auto trim(const std::string_view value) -> std::string_view {
    const auto first = value.find_first_not_of(" \t\r");

    if (first == std::string_view::npos) {
        return {};
    }

    return value.substr(first, value.find_last_not_of(" \t\r") - first + 1);
}

auto element(const std::string_view symbol) -> const core::Element& {
    const auto& table = core::periodic_table();
    if (!table.contains(symbol)) {
        throw ParseError{
            ParseErrorKind::malformed, "unknown element symbol '" + std::string{symbol} + "'", {}};
    }
    return table.element(symbol);
}

auto parse_int(const std::string_view value, const std::string_view field) -> int {
    int result = 0;
    const auto [end, error] = std::from_chars(value.data(), value.data() + value.size(), result);

    if (error != std::errc{} || end != value.data() + value.size()) {
        throw ParseError{ParseErrorKind::malformed, "invalid " + std::string{field}, {}};
    }

    return result;
}

auto parse_trimmed_int(const std::string_view value, const std::string_view field) -> int {
    return parse_int(trim(value), field);
}

auto parse_double(const std::string_view value, const std::string_view field) -> double {
    const auto token = trim(value);

    if (token.empty()) {
        throw ParseError{ParseErrorKind::malformed, "invalid " + std::string{field}, {}};
    }

    // Intentionally not std::from_chars: Apple's libc++ provides it for floating-point values only
    // from macOS 26 (see "Platforms" in DEVELOPMENT.md). Gemmi's locale-independent equivalent also
    // accepts one leading '+'.
    double result = 0.0;
    const auto [end, error] =
        gemmi::fast_from_chars(token.data(), token.data() + token.size(), result);

    if (error != std::errc{} || end != token.data() + token.size()) {
        throw ParseError{ParseErrorKind::malformed, "invalid " + std::string{field}, {}};
    }

    return result;
}

auto fixed_field(const std::string_view line, const std::size_t offset, const std::size_t width,
                 const std::string_view field) -> std::string_view {
    if (line.size() < offset + width) {
        throw ParseError{ParseErrorKind::malformed, "missing " + std::string{field}, {}};
    }

    return line.substr(offset, width);
}

auto bond_order(const std::string_view value, const ::chargefw::io::native::BondFormat format)
    -> core::BondOrder {
    // SYBYL aromatic and amide bond types carry connectivity without an integral order.
    if (format == ::chargefw::io::native::BondFormat::mol2 && (value == "ar" || value == "am")) {
        return core::BondOrder::SINGLE;
    }
    if (format == ::chargefw::io::native::BondFormat::mol) {
        return mol_bond_order(parse_int(value, "MOL bond order"));
    }
    // Dummy, unknown, and not-connected bonds.
    if (value == "du" || value == "un" || value == "nc") {
        throw ParseError{ParseErrorKind::unsupported,
                         "unsupported MOL2 bond type '" + std::string{value} + "'",
                         {}};
    }
    const auto type = parse_int(value, "MOL2 bond type");
    if (type < 1 || type > 4) {
        throw ParseError{
            ParseErrorKind::malformed, "invalid MOL2 bond type '" + std::string{value} + "'", {}};
    }
    return core::bond_order_from_value(type == 4 ? 1 : type);
}

auto mol_bond_order(const int type) -> core::BondOrder {
    // Types 5 to 8 are query bonds, such as "single or double".
    if (type >= 5 && type <= 8) {
        throw ParseError{
            ParseErrorKind::unsupported, "unsupported query bond type " + std::to_string(type), {}};
    }
    if (type < 1 || type > 4) {
        throw ParseError{
            ParseErrorKind::malformed, "invalid MOL bond order " + std::to_string(type), {}};
    }
    return core::bond_order_from_value(type == 4 ? 1 : type);
}

} // namespace chargefw::io::native::common_input
