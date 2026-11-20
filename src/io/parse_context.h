#pragma once

#include <chargefw/io/parse_error.h>

#include <exception>
#include <istream>
#include <new>
#include <utility>

namespace chargefw::io::detail {

// Called from a reader's catch block. A ParseError gets the fields of the location it does not
// have yet: the throw site knows best where the offending data is. A value a molecule cannot hold,
// which the molecule's types refuse with std::invalid_argument, becomes a malformed ParseError.
// Other exceptions, such as a failed read or a broken internal check, pass unchanged.
[[noreturn]] auto rethrow_located(const ParseLocation& location) -> void;

// Throws std::runtime_error for a stream that failed other than by reaching its end, such as a
// file that was not opened. It is not a ParseError: the input was not refused.
auto require_readable(std::istream& input) -> void;

// The result of a call into the parsing library, whose every failure but running out of memory
// means text it cannot parse. A ParseError thrown by ChargeFW's own code in the call passes
// unchanged.
template <typename Parse> [[nodiscard]] auto parsed(Parse&& parse) -> decltype(auto) {
    try {
        return std::forward<Parse>(parse)();
    } catch (const std::bad_alloc&) {
        throw;
    } catch (const ParseError&) {
        throw; // A refusal of ChargeFW's own, which keeps its kind.
    } catch (const std::exception& error) {
        throw ParseError{ParseErrorKind::malformed, error.what(), {}};
    }
}

} // namespace chargefw::io::detail
