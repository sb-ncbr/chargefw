#include <chargefw/adapters/native/pqr_output.h>

#include "adapters/component.h"

#include <chargefw/core/periodic_table.h>

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

namespace chargefw::adapters::native::pqr_output {
namespace {

[[nodiscard]] auto safe_line(std::string value) -> std::string {
    std::ranges::replace_if(
        value, [](const char character) { return character == '\r' || character == '\n'; }, '_');
    return value;
}

// PQR follows PDB, which writes author labels; label_* values are the fallback.
[[nodiscard]] auto author_or_label(const std::optional<std::string>& author,
                                   const std::optional<std::string>& label,
                                   const std::string_view fallback) -> std::string {
    return std::string{detail::first_known_label(author, label, fallback)};
}

[[nodiscard]] auto atom_description(const std::size_t atom_index) -> std::string {
    return "atom " + std::to_string(atom_index + 1);
}

// APBS reads PQR as whitespace-separated words and treats '#' and '%' as comment characters, so a
// label must be one word without them.
auto require_word(const std::string_view value, const std::string_view field,
                  const std::size_t atom_index) -> void {
    if (std::ranges::any_of(value, [](const unsigned char character) {
            return std::isspace(character) != 0 || std::iscntrl(character) != 0 ||
                   character == '#' || character == '%';
        })) {
        throw std::invalid_argument{"PQR output cannot write " + std::string{field} + " '" +
                                    std::string{value} + "' of " + atom_description(atom_index) +
                                    ": it contains whitespace, a control character, '#', or '%'"};
    }
}

// APBS reads the chain only as a single non-numeric character before the residue number; any other
// chain ID shifts the following fields.
auto require_chain(const std::string_view chain, const std::size_t atom_index) -> void {
    require_word(chain, "chain ID", atom_index);
    if (chain.size() != 1 || std::isdigit(static_cast<unsigned char>(chain.front())) != 0) {
        throw std::invalid_argument{"PQR output cannot write chain ID '" + std::string{chain} +
                                    "' of " + atom_description(atom_index) +
                                    ": PQR readers such as APBS require a single non-numeric "
                                    "character"};
    }
}

// PDB files with more than 9,999 residues in a chain write hybrid-36 numbers such as A000, which
// PQR readers would take for a chain ID followed by a number.
auto require_residue_number(const std::string& sequence, const std::size_t atom_index) -> void {
    auto value = 0L;
    const auto* end = sequence.data() + sequence.size();
    const auto [parsed, error] = std::from_chars(sequence.data(), end, value);
    if (error != std::errc{} || parsed != end) {
        throw std::invalid_argument{
            "PQR output cannot write residue number '" + sequence + "' of " +
            atom_description(atom_index) +
            ": it is not an integer; read the structure from mmCIF, which numbers residues in "
            "decimal"};
    }
}

auto require_insertion_code(const std::string_view code, const std::size_t atom_index) -> void {
    if (!code.empty() &&
        (code.size() != 1 || std::isalpha(static_cast<unsigned char>(code.front())) == 0)) {
        throw std::invalid_argument{"PQR output cannot write insertion code '" + std::string{code} +
                                    "' of " + atom_description(atom_index) +
                                    ": it is not a single letter"};
    }
}

auto validate_result(const ChargeCalculationResult& result) -> void {
    if (!result.execution().calculated() || !result.execution().charges.has_value()) {
        throw std::invalid_argument{"PQR output requires a successful calculation"};
    }
    if (result.inputs().size() != 1) {
        throw std::invalid_argument{"PQR output requires exactly one molecule, but the input has " +
                                    std::to_string(result.inputs().size())};
    }

    const auto& record = result.inputs().front();
    if (!record.import_metadata.has_value() ||
        (record.import_metadata->format != MolecularSourceFormat::pdb &&
         record.import_metadata->format != MolecularSourceFormat::mmcif)) {
        throw std::invalid_argument{"PQR output requires a structure read from PDB or mmCIF"};
    }

    const auto& molecule = record.molecule;
    if (molecule.atom_count() == 0) {
        throw std::invalid_argument{"PQR output requires at least one atom"};
    }
    if (molecule.conformer_count() == 0) {
        throw std::invalid_argument{"PQR output requires coordinates"};
    }
    // PQR readers such as APBS ignore MODEL records and would combine all models into one system.
    if (molecule.conformer_count() > 1) {
        throw std::invalid_argument{"PQR output requires one conformer, but the molecule has " +
                                    std::to_string(molecule.conformer_count()) +
                                    "; read the structure with only the first conformer or model"};
    }

    const auto assignments = result.execution().charges->assignments();
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
        throw std::invalid_argument{"PQR output has no van der Waals radius for element " +
                                    std::string{symbol} + " (" + atom_description(atom_index) +
                                    ")"};
    }
    const auto& position = record.molecule.conformer(0)[atom_index];
    if (!std::isfinite(position.x) || !std::isfinite(position.y) || !std::isfinite(position.z)) {
        throw std::invalid_argument{"PQR coordinates must be finite"};
    }
    const auto charge = assignment.charges[atom_index];
    if (!std::isfinite(charge)) {
        throw std::invalid_argument{"PQR charges must be finite"};
    }

    const auto generated_name = std::string{symbol} + std::to_string(atom_index + 1);
    const auto& labels = record.import_metadata->conformers[0].sites[atom_index].structural_labels;
    auto name = labels.has_value()
                    ? author_or_label(labels->author.atom, labels->label.atom, generated_name)
                    : generated_name;
    const auto residue = labels.has_value()
                             ? author_or_label(labels->author.residue, labels->label.residue, "UNL")
                             : std::string{"UNL"};
    const auto chain = labels.has_value()
                           ? author_or_label(labels->author.chain, labels->label.chain, "A")
                           : std::string{"A"};
    const auto sequence =
        labels.has_value() ? author_or_label(labels->author.sequence, labels->label.sequence, "1")
                           : std::string{"1"};
    const auto insertion_code =
        labels.has_value() ? std::string{detail::known_label(labels->insertion_code).value_or("")}
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

auto PqrWriter::write(const ChargeCalculationResult& result,
                      const std::string_view generator_version) const -> void {
    validate_result(result);
    const auto& charge_set = *result.execution().charges;
    const auto& record = result.inputs().front();
    const auto& assignment = charge_set.assignments().front();
    // Every atom is checked before any output, so a refused structure leaves the stream unchanged.
    auto atoms = std::string{};
    for (std::size_t atom_index = 0; atom_index < record.molecule.atom_count(); ++atom_index) {
        append_atom_line(atoms, record, assignment, atom_index);
    }

    std::print(*output_, "REMARK   1 Generated by ChargeFW");
    if (!generator_version.empty()) {
        std::print(*output_, " {}", safe_line(std::string{generator_version}));
    }
    std::print(*output_, "; method={}", safe_line(std::string{charge_set.method_id()}));
    if (const auto parameter_set_id = charge_set.parameter_set_id(); parameter_set_id.has_value()) {
        std::print(*output_, "; parameter_set={}", safe_line(std::string{*parameter_set_id}));
    }
    std::println(*output_);
    std::println(*output_, "REMARK   1 Radii are van der Waals radii of the elements");
    std::print(*output_, "{}", atoms);
    std::println(*output_, "END");
    if (!*output_) {
        throw std::runtime_error{"failed to write PQR output"};
    }
}

} // namespace chargefw::adapters::native::pqr_output
