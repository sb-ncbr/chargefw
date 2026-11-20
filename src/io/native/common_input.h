#pragma once

#include <chargefw/core/bond.h>
#include <chargefw/core/periodic_table.h>
#include <chargefw/io/molecule_record.h>
#include <chargefw/io/parse_error.h>

#include "bond_format.h"

#include <cstddef>
#include <istream>
#include <optional>
#include <string>
#include <string_view>

namespace chargefw::io::native::common_input {

// Called from a parser's catch block. A ParseError gets the record and line it does not have yet;
// a throw site inside a record leaves its location empty unless the offending data is elsewhere.
// A value a molecule cannot hold, which the molecule's types refuse with std::invalid_argument,
// becomes a malformed ParseError. Other exceptions, such as a failed read, pass unchanged.
[[noreturn]] auto rethrow_with_context(const MoleculeRecordIdentity& identity,
                                       std::optional<std::size_t> line) -> void;

// Whether the input has no more characters; a stream that failed otherwise throws
// std::runtime_error, which is not a ParseError: the input was not refused.
[[nodiscard]] auto at_end(std::istream& input) -> bool;

// The next line, or nothing at the end of the input; a failed stream throws as in at_end().
[[nodiscard]] auto try_read_line(std::istream& input, std::size_t& line)
    -> std::optional<std::string>;

[[nodiscard]] auto read_line(std::istream& input, std::size_t& line, std::string_view record_name)
    -> std::string;

[[nodiscard]] auto trim(std::string_view value) -> std::string_view;

[[nodiscard]] auto parse_int(std::string_view value, std::string_view field) -> int;

// The element with this symbol; an unknown symbol is input that cannot be read, not a bad index.
[[nodiscard]] auto element(std::string_view symbol) -> const core::Element&;

[[nodiscard]] auto parse_trimmed_int(std::string_view value, std::string_view field) -> int;

[[nodiscard]] auto parse_double(std::string_view value, std::string_view field) -> double;

[[nodiscard]] auto fixed_field(std::string_view line, std::size_t offset, std::size_t width,
                               std::string_view field) -> std::string_view;

// The order of a bond type of the format; aromatic bonds are single. A type the format defines
// but a molecule cannot hold, such as a query bond, is unsupported; any other is malformed.
[[nodiscard]] auto bond_order(std::string_view value, ::chargefw::io::native::BondFormat format)
    -> core::BondOrder;
[[nodiscard]] auto mol_bond_order(int type) -> core::BondOrder;

} // namespace chargefw::io::native::common_input
