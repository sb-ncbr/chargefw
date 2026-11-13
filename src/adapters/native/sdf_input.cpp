#include <chargefw/adapters/native/sdf_input.h>

#include "common_input.h"
#include "mol_record_parser.h"

#include <istream>
#include <string>
#include <utility>

namespace chargefw::adapters::native::sdf_input {
namespace {

// Blank lines may begin a record with an empty title, so whitespace is skipped only to detect
// trailing whitespace at the end of the file and is otherwise restored.
auto at_end_of_records(std::istream& input) -> bool {
    if (input.eof()) {
        return true;
    }
    const auto start = input.tellg();
    if (start == std::istream::pos_type(-1)) {
        return input.peek() == std::char_traits<char>::eof();
    }

    input >> std::ws;
    if (input.eof()) {
        return true;
    }
    input.seekg(start);
    return false;
}

auto consume_to_sdf_delimiter(std::istream& input, std::size_t& line_number) -> void {
    while (const auto line = common_input::try_read_line(input, line_number)) {
        if (*line == "$$$$") {
            return;
        }
    }
}

} // namespace

SdfReader::SdfReader(std::istream& input, std::string source_name)
    : input_{std::addressof(input)}, source_{std::move(source_name)} {}

auto SdfReader::next() -> std::optional<ImportedMoleculeRecord> {
    if (at_end_of_records(*input_)) {
        return std::nullopt;
    }

    const auto identity = MoleculeRecordIdentity{
        .source_name = source_, .record_index = record_index_++, .record_id = {}};
    auto result = mol_input::parse_mol(*input_, identity, line_);
    result.import_metadata->format = MolecularSourceFormat::sdf;
    try {
        consume_to_sdf_delimiter(*input_, line_);
    } catch (...) {
        common_input::rethrow_with_context(result.identity, line_);
    }
    return result;
}

} // namespace chargefw::adapters::native::sdf_input
