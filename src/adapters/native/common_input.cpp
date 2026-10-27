#include "common_input.h"

#include <chargefw/core/molecule.h>

#include <gemmi/atof.hpp>

#include <charconv>
#include <stdexcept>
#include <string>
#include <utility>

namespace chargefw::adapters::native::common_input {

auto read_line(std::istream& input, std::size_t& line, const std::string_view record_name)
    -> std::string {
    std::string result;

    if (!std::getline(input, result)) {
        throw std::runtime_error{"unexpected end of " + std::string{record_name} + " record"};
    }

    if (result.ends_with('\r')) {
        result.pop_back();
    }

    ++line;
    return result;
}

auto trim(const std::string_view value) -> std::string_view {
    const auto first = value.find_first_not_of(" \t\r");

    if (first == std::string_view::npos) {
        return {};
    }

    return value.substr(first, value.find_last_not_of(" \t\r") - first + 1);
}

auto parse_int(const std::string_view value, const std::string_view field) -> int {
    int result = 0;
    const auto [end, error] = std::from_chars(value.data(), value.data() + value.size(), result);

    if (error != std::errc{} || end != value.data() + value.size()) {
        throw std::runtime_error{"invalid " + std::string{field}};
    }

    return result;
}

auto parse_trimmed_int(const std::string_view value, const std::string_view field) -> int {
    return parse_int(trim(value), field);
}

auto parse_double(const std::string_view value, const std::string_view field) -> double {
    const auto token = trim(value);

    if (token.empty()) {
        throw std::runtime_error{"invalid " + std::string{field}};
    }

    // Intentionally not std::from_chars: Apple's libc++ provides it for floating-point values only
    // from macOS 26 (see "Platforms" in DEVELOPMENT.md). Gemmi's locale-independent equivalent also
    // accepts one leading '+'.
    double result = 0.0;
    const auto [end, error] =
        gemmi::fast_from_chars(token.data(), token.data() + token.size(), result);

    if (error != std::errc{} || end != token.data() + token.size()) {
        throw std::runtime_error{"invalid " + std::string{field}};
    }

    return result;
}

auto parse_trimmed_double(const std::string_view value, const std::string_view field) -> double {
    return parse_double(trim(value), field);
}

auto fixed_field(const std::string_view line, const std::size_t offset, const std::size_t width,
                 const std::string_view field) -> std::string_view {
    if (line.size() < offset + width) {
        throw std::runtime_error{"missing " + std::string{field}};
    }

    return line.substr(offset, width);
}

auto bond_order(const std::string_view value, const ::chargefw::adapters::native::BondFormat format)
    -> core::BondOrder {
    // SYBYL aromatic and amide bond types carry connectivity without an integral order.
    if (format == ::chargefw::adapters::native::BondFormat::mol2 &&
        (value == "ar" || value == "am")) {
        return core::BondOrder::SINGLE;
    }
    const auto numeric_value = parse_int(
        value, format == ::chargefw::adapters::native::BondFormat::mol ? "MOL bond order"
                                                                       : "MOL2 bond type");
    return core::bond_order_from_value(numeric_value == 4 ? 1 : numeric_value);
}

auto make_record(std::vector<core::Atom> atoms, std::vector<core::Bond> bonds,
                 std::vector<core::Conformer> conformers, MoleculeRecordIdentity identity,
                 std::string name, std::vector<MoleculeRecordDiagnostic> diagnostics,
                 std::optional<MoleculeImportMetadata> import_metadata) -> ImportedMoleculeRecord {
    if (name.empty()) {
        name = identity.record_id.display_string();
    }

    return ImportedMoleculeRecord{.molecule =
                                      core::Molecule{std::move(atoms), std::move(bonds),
                                                     std::move(conformers), std::move(name)},
                                  .identity = std::move(identity),
                                  .caller_atom_ids = std::nullopt,
                                  .diagnostics = std::move(diagnostics),
                                  .import_metadata = std::move(import_metadata)};
}

} // namespace chargefw::adapters::native::common_input
