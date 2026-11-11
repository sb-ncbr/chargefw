#include <chargefw/adapters/native/mol2_input.h>

#include "common_input.h"

#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/periodic_table.h>
#include <chargefw/core/position.h>

#include <cstddef>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <vector>

namespace chargefw::adapters::native::mol2_input {
namespace {

namespace common = chargefw::adapters::native::common_input;

constexpr std::string_view molecule_marker{"@<TRIPOS>MOLECULE"};
constexpr std::string_view atom_marker{"@<TRIPOS>ATOM"};
constexpr std::string_view bond_marker{"@<TRIPOS>BOND"};
constexpr std::string_view unity_atom_attributes_marker{"@<TRIPOS>UNITY_ATOM_ATTR"};

[[nodiscard]] auto element_symbol(const std::string_view atom_type) -> std::string_view {
    const auto separator = atom_type.find('.');
    const auto symbol = atom_type.substr(0, separator);

    if (symbol.empty() || symbol == "Du" || symbol == "LP" || symbol == "Any" || symbol == "Hal" ||
        symbol == "Het" || symbol == "Hev") {
        throw std::runtime_error{"unsupported MOL2 atom type '" + std::string{atom_type} + "'"};
    }

    return symbol;
}

// Skips the data lines of a section and returns the next record type indicator, or nothing at the
// end of the input.
[[nodiscard]] auto next_section(std::istream& input, std::size_t& line)
    -> std::optional<std::string> {
    while (auto text = common::try_read_line(input, line)) {
        if (text->starts_with('@')) {
            return text;
        }
    }
    return std::nullopt;
}

// The atoms of an ATOM section, all neutral; formal charges are applied from UNITY_ATOM_ATTR.
struct AtomSection {
    std::vector<core::Atom> atoms;
    std::vector<core::Position> positions;
    std::vector<SourceAtomReference> references;
    std::unordered_map<int, std::size_t> indices;
    std::optional<std::size_t> partial_charge_line;
};

[[nodiscard]] auto read_atoms(std::istream& input, std::size_t& line, const int atom_count)
    -> AtomSection {
    auto result = AtomSection{};

    for (int index = 0; index < atom_count; ++index) {
        std::istringstream atom_line{common::read_line(input, line, "MOL2")};
        std::string source_id_text;
        std::string atom_name;
        std::string x_text;
        std::string y_text;
        std::string z_text;
        std::string atom_type;
        std::string subst_id;
        std::string subst_name;
        std::string charge_text;

        if (!(atom_line >> source_id_text >> atom_name >> x_text >> y_text >> z_text >>
              atom_type)) {
            throw std::runtime_error{"invalid MOL2 atom record"};
        }

        static_cast<void>(atom_line >> subst_id >> subst_name >> charge_text);
        const auto source_id = common::parse_int(source_id_text, "MOL2 atom id");

        if (source_id <= 0 || result.indices.contains(source_id)) {
            throw std::runtime_error{"duplicate or invalid MOL2 atom id"};
        }

        if (!charge_text.empty()) {
            const auto charge = common::parse_double(charge_text, "MOL2 atom charge");

            if (charge != 0.0 && !result.partial_charge_line.has_value()) {
                result.partial_charge_line = line;
            }
        }

        const auto& element = core::periodic_table().element(element_symbol(atom_type));
        result.indices.emplace(source_id, result.atoms.size());
        result.atoms.emplace_back(element.atomic_number, 0, std::move(atom_name));
        result.references.push_back(SourceAtomReference{.position = static_cast<std::size_t>(index),
                                                        .id = std::move(source_id_text),
                                                        .structural_labels = nullptr});
        result.positions.push_back(
            core::Position{.x = common::parse_double(x_text, "MOL2 x coordinate"),
                           .y = common::parse_double(y_text, "MOL2 y coordinate"),
                           .z = common::parse_double(z_text, "MOL2 z coordinate")});
    }

    return result;
}

// A bond of a BOND section, by source atom IDs.
struct SourceBond {
    int first;
    int second;
    std::string type;
    std::size_t line;
};

struct FormalCharge {
    int value;
    std::size_t line;
};

[[nodiscard]] auto read_bonds(std::istream& input, std::size_t& line, const int bond_count)
    -> std::vector<SourceBond> {
    auto result = std::vector<SourceBond>{};

    for (int index = 0; index < bond_count; ++index) {
        std::istringstream bond_line{common::read_line(input, line, "MOL2")};
        std::string bond_id;
        std::string first_text;
        std::string second_text;
        std::string type;

        if (!(bond_line >> bond_id >> first_text >> second_text >> type)) {
            throw std::runtime_error{"invalid MOL2 bond record"};
        }

        result.push_back(SourceBond{.first = common::parse_int(first_text, "MOL2 bond atom"),
                                    .second = common::parse_int(second_text, "MOL2 bond atom"),
                                    .type = std::move(type),
                                    .line = line});
    }

    return result;
}

// Reads the formal charges of a UNITY_ATOM_ATTR section by atom ID and returns the next record type
// indicator. Each record is an atom ID and an attribute count, followed by one line per attribute
// with its name and an optional value; "charge" is the formal charge.
[[nodiscard]] auto read_formal_charges(std::istream& input, std::size_t& line,
                                       std::unordered_map<int, FormalCharge>& charges)
    -> std::optional<std::string> {
    while (auto text = common::try_read_line(input, line)) {
        if (text->starts_with('@')) {
            return text;
        }
        std::istringstream header{*text};
        std::string atom_id_text;
        std::string attribute_count_text;
        if (text->starts_with('#') || !(header >> atom_id_text)) {
            continue;
        }
        if (!(header >> attribute_count_text)) {
            throw std::runtime_error{"invalid MOL2 UNITY_ATOM_ATTR record"};
        }
        const auto atom_id = common::parse_int(atom_id_text, "MOL2 UNITY_ATOM_ATTR atom");
        const auto atom_line = line;
        const auto attribute_count =
            common::parse_int(attribute_count_text, "MOL2 UNITY_ATOM_ATTR count");
        if (attribute_count < 0) {
            throw std::runtime_error{"invalid MOL2 UNITY_ATOM_ATTR record"};
        }
        for (int attribute = 0; attribute < attribute_count; ++attribute) {
            const auto attribute_text = common::read_line(input, line, "MOL2");
            std::istringstream attribute_line{attribute_text};
            std::string attribute_name;
            std::string value;
            if (attribute_text.starts_with('@') || !(attribute_line >> attribute_name)) {
                throw std::runtime_error{"MOL2 UNITY_ATOM_ATTR record has fewer attributes "
                                         "than its count"};
            }
            attribute_line >> value;
            if (attribute_name == "charge") {
                charges[atom_id] = FormalCharge{
                    common::parse_int(value, "MOL2 UNITY_ATOM_ATTR charge"), atom_line};
            }
        }
    }
    return std::nullopt;
}

// Reads one record. Its sections may come in any order; ChargeFW uses ATOM, BOND, and the formal
// charges in UNITY_ATOM_ATTR, and skips the others. Sets next_record when the record ended at the
// MOLECULE line of the next one.
[[nodiscard]] auto parse_record(std::istream& input, MoleculeRecordIdentity identity,
                                bool& next_record, std::size_t& line) -> ImportedMoleculeRecord {
    std::optional<std::size_t> validation_line;
    try {
        const auto name = common::read_line(input, line, "MOL2");
        if (identity.record_id.empty()) {
            identity.record_id = name;
        }
        const auto counts_line = common::read_line(input, line, "MOL2");
        std::istringstream counts{counts_line};
        std::string atom_count_text;
        std::string bond_count_text;

        if (!(counts >> atom_count_text >> bond_count_text)) {
            throw std::runtime_error{"invalid MOL2 counts"};
        }

        const auto atom_count = common::parse_int(atom_count_text, "MOL2 atom count");
        const auto bond_count = common::parse_int(bond_count_text, "MOL2 bond count");

        if (atom_count <= 0 || bond_count < 0) {
            throw std::runtime_error{"invalid MOL2 counts"};
        }

        auto atom_section = std::optional<AtomSection>{};
        auto source_bonds = std::optional<std::vector<SourceBond>>{};
        auto formal_charges = std::unordered_map<int, FormalCharge>{};
        auto section = next_section(input, line);
        while (section.has_value() && common::trim(*section) != molecule_marker) {
            const auto marker = common::trim(*section);
            if (marker == atom_marker) {
                if (atom_section.has_value()) {
                    throw std::runtime_error{"repeated MOL2 ATOM section"};
                }
                atom_section = read_atoms(input, line, atom_count);
                section = next_section(input, line);
            } else if (marker == bond_marker) {
                if (source_bonds.has_value()) {
                    throw std::runtime_error{"repeated MOL2 BOND section"};
                }
                source_bonds = read_bonds(input, line, bond_count);
                section = next_section(input, line);
            } else if (marker == unity_atom_attributes_marker) {
                section = read_formal_charges(input, line, formal_charges);
            } else {
                section = next_section(input, line);
            }
        }
        next_record = section.has_value();

        if (!atom_section.has_value()) {
            throw std::runtime_error{"missing MOL2 ATOM section"};
        }
        if (!source_bonds.has_value()) {
            throw std::runtime_error{"expected MOL2 BOND section"};
        }
        auto& [atoms, positions, atom_references, atom_indices, partial_charge_line] =
            *atom_section;

        for (const auto& [atom_id, charge] : formal_charges) {
            validation_line = charge.line;
            const auto index = atom_indices.find(atom_id);
            if (index == atom_indices.end()) {
                throw std::runtime_error{"MOL2 UNITY_ATOM_ATTR references an unknown atom"};
            }
            auto& atom = atoms[index->second];
            atom = core::Atom{atom.atomic_number(), charge.value, std::string{atom.name()}};
        }

        std::vector<core::Bond> bonds;
        bonds.reserve(source_bonds->size());
        auto aromatic_bonds = std::vector<std::size_t>{}; // Indices in bonds.

        for (const auto& [first_id, second_id, type, source_line] : *source_bonds) {
            validation_line = source_line;
            const auto first = atom_indices.find(first_id);
            const auto second = atom_indices.find(second_id);

            if (first == atom_indices.end() || second == atom_indices.end()) {
                throw std::runtime_error{"MOL2 bond references an unknown atom"};
            }

            bonds.emplace_back(first->second, second->second,
                               common::bond_order(type, BondFormat::mol2));
            if (type == "ar" || type == "4") {
                aromatic_bonds.push_back(bonds.size() - 1);
            }
        }

        validation_line.reset();

        auto issues = std::vector<ImportIssue>{};

        if (partial_charge_line.has_value()) {
            issues.push_back(ImportIssue{
                .kind = ImportIssueKind::partial_charges_ignored,
                .message = "MOL2 partial charges were ignored; they are not formal charges and "
                           "were not used for calculation.",
                .line = partial_charge_line});
        }
        if (!aromatic_bonds.empty()) {
            issues.push_back(common::aromatic_bonds_issue(std::move(aromatic_bonds)));
        }

        auto sites = atom_references;
        auto metadata = MoleculeImportMetadata{
            .format = MolecularSourceFormat::mol2,
            .atoms = std::move(atom_references),
            .conformers = {SourceConformerReference{
                .position = 0, .id = std::nullopt, .sites = std::move(sites)}},
            .record_selection = std::nullopt,
            .alternate_location_selection = std::nullopt,
            .conformer_selection = "all",
            .bond_strategy = std::nullopt,
            .source_connectivity = bond_count == 0 ? SourceConnectivity::explicitly_empty
                                                   : SourceConnectivity::present,
            .components = {},
        };
        return common::make_record(std::move(atoms), std::move(bonds),
                                   {core::Conformer{std::move(positions), "input"}}, identity, {},
                                   std::move(issues), std::move(metadata));
    } catch (...) {
        common::rethrow_with_context(identity, validation_line.value_or(line));
    }
}

} // namespace

Mol2Reader::Mol2Reader(std::istream& input, std::string source)
    : input_{std::addressof(input)}, source_{std::move(source)} {}

auto Mol2Reader::next() -> std::optional<ImportedMoleculeRecord> {
    auto found = std::exchange(next_record_started_, false);
    while (!found) {
        const auto line = common::try_read_line(*input_, line_);
        if (!line.has_value()) {
            return std::nullopt;
        }
        found = line->starts_with('@') && common::trim(*line) == molecule_marker;
    }
    return parse_record(
        *input_,
        MoleculeRecordIdentity{.source = source_, .record_index = record_index_++, .record_id = {}},
        next_record_started_, line_);
}

} // namespace chargefw::adapters::native::mol2_input
