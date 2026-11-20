#include <chargefw/io/pqr_output.h>

#include "io/component.h"
#include "io/output_helpers.h"

#include <chargefw/core/periodic_table.h>
#include <chargefw/io/output_error.h>

#include <algorithm>
#include <cctype>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <format>
#include <iterator>
#include <memory>
#include <optional>
#include <ostream>
#include <print>
#include <stdexcept>
#include <string>
#include <string_view>

namespace chargefw::io {
namespace {

// PQR follows PDB, which writes author labels; label_* values are the fallback.
[[nodiscard]] auto author_or_label(const std::optional<std::string>& author,
                                   const std::optional<std::string>& label,
                                   const std::string_view fallback) -> std::string {
    return std::string{detail::first_known_label(author, label, fallback)};
}

// Where a refusal of one atom is. PQR output has one molecule, at index 0.
[[nodiscard]] auto at_atom(const std::size_t atom_index) -> core::IssueLocation {
    return core::IssueLocation{.atom_index = atom_index, .atom_indices = {atom_index}};
}

// APBS reads PQR as whitespace-separated words and treats '#' and '%' as comment characters, so a
// label must be one word without them.
auto require_word(const std::string_view value, const std::string_view field,
                  const std::size_t atom_index) -> void {
    if (std::ranges::any_of(value, [](const unsigned char character) {
            return std::isspace(character) != 0 || std::iscntrl(character) != 0 ||
                   character == '#' || character == '%';
        })) {
        throw OutputError{OutputErrorKind::unrepresentable,
                          "PQR output cannot write " + std::string{field} + " '" +
                              std::string{value} +
                              "': it contains whitespace, a control character, '#', or '%'",
                          0, at_atom(atom_index)};
    }
}

// APBS reads the chain only as a single non-numeric character before the residue number; any other
// chain ID shifts the following fields.
auto require_chain(const std::string_view chain, const std::size_t atom_index) -> void {
    require_word(chain, "chain ID", atom_index);
    if (chain.size() != 1 || std::isdigit(static_cast<unsigned char>(chain.front())) != 0) {
        throw OutputError{OutputErrorKind::unrepresentable,
                          "PQR output cannot write chain ID '" + std::string{chain} +
                              "': PQR readers such as APBS require a single non-numeric character",
                          0, at_atom(atom_index)};
    }
}

// PDB files with more than 9,999 residues in a chain write hybrid-36 numbers such as A000, which
// PQR readers would take for a chain ID followed by a number.
auto require_residue_number(const std::string& sequence, const std::size_t atom_index) -> void {
    auto value = 0L;
    const auto* end = sequence.data() + sequence.size();
    const auto [parsed, error] = std::from_chars(sequence.data(), end, value);
    if (error != std::errc{} || parsed != end) {
        throw OutputError{
            OutputErrorKind::unrepresentable,
            "PQR output cannot write residue number '" + sequence +
                "': it is not an integer; read the structure from mmCIF, which numbers residues in "
                "decimal",
            0, at_atom(atom_index)};
    }
}

auto require_insertion_code(const std::string_view code, const std::size_t atom_index) -> void {
    if (!code.empty() &&
        (code.size() != 1 || std::isalpha(static_cast<unsigned char>(code.front())) == 0)) {
        throw OutputError{OutputErrorKind::unrepresentable,
                          "PQR output cannot write insertion code '" + std::string{code} +
                              "': it is not a single letter",
                          0, at_atom(atom_index)};
    }
}

auto validate_result(const ChargeCalculationResult& result) -> void {
    if (result.charges() == nullptr) {
        throw std::invalid_argument{"PQR output requires a successful calculation"};
    }
    // The source first: input that is not a structure is never written as PQR, however many
    // molecules it has.
    if (!std::ranges::all_of(result.inputs(), [](const ImportedMoleculeRecord& record) {
            return record.import_metadata.has_value() &&
                   (record.import_metadata->format == MolecularSourceFormat::pdb ||
                    record.import_metadata->format == MolecularSourceFormat::mmcif);
        })) {
        throw std::invalid_argument{"PQR output requires a structure read from PDB or mmCIF"};
    }
    if (result.inputs().size() != 1) {
        throw OutputError{OutputErrorKind::unsupported_shape,
                          "PQR output requires exactly one molecule, but the input has " +
                              std::to_string(result.inputs().size())};
    }

    const auto& record = result.inputs().front();

    const auto& molecule = record.molecule;
    if (molecule.atom_count() == 0) {
        throw OutputError{OutputErrorKind::unsupported_shape,
                          "PQR output requires at least one atom", 0};
    }
    if (molecule.conformer_count() == 0) {
        throw OutputError{OutputErrorKind::missing_coordinates, "PQR output requires coordinates",
                          0};
    }
    // PQR readers such as APBS ignore MODEL records and would combine all models into one system.
    if (molecule.conformer_count() > 1) {
        throw OutputError{OutputErrorKind::unsupported_shape,
                          "PQR output requires one conformer, but the molecule has " +
                              std::to_string(molecule.conformer_count()) +
                              "; read the structure with only the first conformer or model",
                          0};
    }

    const auto assignments = result.charges()->assignments();
    if (assignments.size() != 1 || assignments.front().target.molecule_index != 0 ||
        assignments.front().charges.size() != molecule.atom_count()) {
        throw std::invalid_argument{"PQR charge assignment does not match the molecule"};
    }
}

auto append_atom_line(std::string& text, const ImportedMoleculeRecord& record,
                      const charges::ChargeAssignment& assignment, const std::size_t atom_index)
    -> void {
    const auto atomic_number = record.molecule.atom(atom_index).atomic_number();
    const auto symbol = core::element_symbol(atomic_number);
    const auto radius = core::periodic_table().element(atomic_number).van_der_waals_radius;
    if (!(radius > 0.0)) {
        throw OutputError{OutputErrorKind::unrepresentable,
                          "PQR output has no van der Waals radius for element " +
                              std::string{symbol},
                          0, at_atom(atom_index)};
    }
    const auto& position = record.molecule.conformer(0)[atom_index];
    if (!std::isfinite(position.x) || !std::isfinite(position.y) || !std::isfinite(position.z)) {
        throw OutputError{OutputErrorKind::unrepresentable, "PQR coordinates must be finite", 0,
                          at_atom(atom_index)};
    }
    const auto charge = assignment.charges[atom_index];
    if (!std::isfinite(charge)) {
        throw OutputError{OutputErrorKind::unrepresentable, "PQR charges must be finite", 0,
                          at_atom(atom_index)};
    }

    const auto generated_name = std::string{symbol} + std::to_string(atom_index + 1);
    const auto* labels =
        record.import_metadata->conformers[0].sites[atom_index].structural_labels.get();
    auto name = labels != nullptr
                    ? author_or_label(labels->author.atom, labels->label.atom, generated_name)
                    : generated_name;
    const auto residue = labels != nullptr
                             ? author_or_label(labels->author.residue, labels->label.residue, "UNL")
                             : std::string{"UNL"};
    const auto chain = labels != nullptr
                           ? author_or_label(labels->author.chain, labels->label.chain, "A")
                           : std::string{"A"};
    const auto sequence =
        labels != nullptr ? author_or_label(labels->author.sequence, labels->label.sequence, "1")
                          : std::string{"1"};
    const auto insertion_code =
        labels != nullptr ? std::string{detail::known_label(labels->insertion_code).value_or("")}
                          : std::string{};
    require_word(name, "atom name", atom_index);
    require_word(residue, "residue name", atom_index);
    require_chain(chain, atom_index);
    require_residue_number(sequence, atom_index);
    require_insertion_code(insertion_code, atom_index);

    // PDB starts names of one-letter elements in column 14 unless they fill four columns.
    if (symbol.size() == 1 && name.size() < 4) {
        name.insert(0, 1, ' ');
    }
    // Every field is separated by a space, unlike PDB columns, so that readers splitting on
    // whitespace, such as APBS, see the same fields for residue numbers above 999.
    std::format_to(std::back_inserter(text),
                   "ATOM  {:>5} {:<4} {:>3} {} {:>4}{:1} {:8.3f} {:8.3f} {:8.3f} {:7.4f} {:6.4f}\n",
                   atom_index + 1, name, residue, chain, sequence, insertion_code, position.x,
                   position.y, position.z, charge, radius);
}

} // namespace

PqrWriter::PqrWriter(std::ostream& output) : output_{std::addressof(output)} {}

auto PqrWriter::write(const ChargeCalculationResult& result) const -> void {
    validate_result(result);
    const auto& charge_set = *result.charges();
    const auto& record = result.inputs().front();
    const auto& assignment = charge_set.assignments().front();
    // Every atom is checked before any output, so a refused structure leaves the stream unchanged.
    auto atoms = std::string{};
    for (std::size_t atom_index = 0; atom_index < record.molecule.atom_count(); ++atom_index) {
        append_atom_line(atoms, record, assignment, atom_index);
    }

    std::print(*output_, "REMARK   1 ");
    detail::write_generator_comment(*output_, charge_set);
    std::println(*output_, "REMARK   1 Radii are van der Waals radii of the elements");
    std::print(*output_, "{}", atoms);
    std::println(*output_, "END");
    if (!*output_) {
        throw std::runtime_error{"failed to write PQR output"};
    }
}

} // namespace chargefw::io
