#include <chargefw/adapters/native/mol_input.h>

#include "common_input.h"

#include "adapters/imported_records.h"
#include "mol_record_parser.h"

#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/position.h>

#include <cstddef>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace chargefw::adapters::native::mol_input {
namespace {

namespace common = chargefw::adapters::native::common_input;

constexpr std::string_view v30_prefix{"M  V30 "};

[[nodiscard]] auto atom_from_symbol(std::string_view symbol, const int formal_charge)
    -> core::Atom {
    const auto normalized = common::trim(symbol);

    if (normalized.empty() || normalized == "*" || normalized == "A" || normalized == "Q" ||
        normalized == "L" || normalized == "LP" || normalized == "R" || normalized == "R#") {
        throw std::runtime_error{"unsupported query atom '" + std::string{normalized} + "'"};
    }

    // Molfiles have no atom names; the symbol is the element, not a name.
    const auto& element = common::element(normalized);
    return core::Atom{element.atomic_number, formal_charge};
}

// The charge field of a V2000 atom line; writers may end the line after the atom symbol.
[[nodiscard]] auto v2000_atom_charge(const std::string_view atom_line) -> int {
    if (atom_line.size() < 39) {
        return 0;
    }

    const auto code = common::trim(atom_line.substr(36, 3));
    if (code.empty()) {
        return 0;
    }

    switch (common::parse_int(code, "V2000 atom charge")) {
    case 0:
    case 4: // doublet radical, uncharged
        return 0;
    case 1:
        return 3;
    case 2:
        return 2;
    case 3:
        return 1;
    case 5:
        return -1;
    case 6:
        return -2;
    case 7:
        return -3;
    default:
        throw std::runtime_error{"invalid V2000 atom charge '" + std::string{code} + "'"};
    }
}

struct ParsedMolecule {
    std::vector<core::Atom> atoms;
    std::vector<core::Bond> bonds;
    std::vector<core::Position> positions;
    std::vector<SourceAtomReference> atom_references;
    SourceConnectivity source_connectivity = SourceConnectivity::absent;
    std::vector<ImportIssue> issues;
    std::vector<std::size_t> aromatic_bonds; // Indices in bonds.
};

auto parse_v2000(std::istream& input, const std::string_view counts, std::size_t& line)
    -> ParsedMolecule {
    const auto atom_count = common::parse_trimmed_int(
        common::fixed_field(counts, 0, 3, "V2000 atom count"), "V2000 atom count");
    const auto bond_count = common::parse_trimmed_int(
        common::fixed_field(counts, 3, 3, "V2000 bond count"), "V2000 bond count");

    if (atom_count <= 0 || bond_count < 0) {
        throw std::runtime_error{"invalid V2000 counts"};
    }

    auto result = ParsedMolecule{};
    result.atoms.reserve(static_cast<std::size_t>(atom_count));
    result.positions.reserve(static_cast<std::size_t>(atom_count));
    result.atom_references.reserve(static_cast<std::size_t>(atom_count));
    result.bonds.reserve(static_cast<std::size_t>(bond_count));
    result.source_connectivity =
        bond_count == 0 ? SourceConnectivity::explicitly_empty : SourceConnectivity::present;
    auto charges_from_properties = false;

    for (int index = 0; index < atom_count; ++index) {
        const auto atom_line = common::read_line(input, line, "MOL");
        const auto x = common::parse_double(
            common::fixed_field(atom_line, 0, 10, "V2000 x coordinate"), "V2000 x coordinate");
        const auto y = common::parse_double(
            common::fixed_field(atom_line, 10, 10, "V2000 y coordinate"), "V2000 y coordinate");
        const auto z = common::parse_double(
            common::fixed_field(atom_line, 20, 10, "V2000 z coordinate"), "V2000 z coordinate");
        const auto symbol = common::fixed_field(atom_line, 31, 3, "V2000 atom symbol");

        result.atoms.push_back(atom_from_symbol(symbol, v2000_atom_charge(atom_line)));
        result.positions.push_back(core::Position{.x = x, .y = y, .z = z});
        result.atom_references.push_back(
            SourceAtomReference{.position = static_cast<std::size_t>(index),
                                .id = std::nullopt,
                                .structural_labels = nullptr});
    }

    for (int index = 0; index < bond_count; ++index) {
        const auto bond_line = common::read_line(input, line, "MOL");
        const auto first = common::parse_trimmed_int(
            common::fixed_field(bond_line, 0, 3, "V2000 bond atom"), "V2000 bond atom");
        const auto second = common::parse_trimmed_int(
            common::fixed_field(bond_line, 3, 3, "V2000 bond atom"), "V2000 bond atom");
        const auto order = common::trim(common::fixed_field(bond_line, 6, 3, "V2000 bond order"));

        if (first <= 0 || second <= 0 || first > atom_count || second > atom_count) {
            throw std::runtime_error{"V2000 bond references an unknown atom"};
        }

        result.bonds.emplace_back(static_cast<std::size_t>(first - 1),
                                  static_cast<std::size_t>(second - 1),
                                  common::bond_order(order, BondFormat::mol));
        if (order == "4") {
            result.aromatic_bonds.push_back(result.bonds.size() - 1);
        }
    }

    // Of the lines after the bonds, only the charge lines are read; the rest, such as isotopes,
    // S-groups, and atom aliases, is ignored up to M  END.
    while (true) {
        const auto property_line = common::read_line(input, line, "MOL");

        if (property_line == "M  END") {
            return result;
        }
        if (property_line == "$$$$") {
            throw std::runtime_error{"missing M  END"};
        }
        // An atom alias (A) or group abbreviation (G) line is followed by its text, which is
        // skipped unread even if it looks like a property line.
        if (property_line.starts_with("A  ") || property_line.starts_with("G  ")) {
            if (common::read_line(input, line, "MOL") == "$$$$") {
                throw std::runtime_error{"missing M  END"};
            }
            continue;
        }

        // M  CHG and M  RAD lines replace all charges and radicals of the atom block.
        if ((property_line.starts_with("M  CHG") || property_line.starts_with("M  RAD")) &&
            !charges_from_properties) {
            charges_from_properties = true;
            for (auto& atom : result.atoms) {
                atom = core::Atom{atom.atomic_number(), 0, std::string{atom.name()}};
            }
        }

        if (!property_line.starts_with("M  CHG")) {
            continue;
        }

        const auto count = common::parse_trimmed_int(
            common::fixed_field(property_line, 6, 3, "M  CHG count"), "M  CHG count");

        if (count < 0 || property_line.size() < 9 + static_cast<std::size_t>(count) * 8) {
            throw std::runtime_error{"invalid M  CHG record"};
        }

        for (int index = 0; index < count; ++index) {
            const auto offset = 9 + static_cast<std::size_t>(index) * 8;
            const auto atom_number = common::parse_trimmed_int(
                common::fixed_field(property_line, offset, 4, "M  CHG atom"), "M  CHG atom");
            const auto charge = common::parse_trimmed_int(
                common::fixed_field(property_line, offset + 4, 4, "M  CHG charge"),
                "M  CHG charge");

            if (atom_number <= 0 || atom_number > atom_count) {
                throw std::runtime_error{"M  CHG references an unknown atom"};
            }

            const auto atom_index = static_cast<std::size_t>(atom_number - 1);
            const auto& atom = result.atoms[atom_index];
            result.atoms[atom_index] =
                core::Atom{atom.atomic_number(), charge, std::string{atom.name()}};
        }
    }
}

[[nodiscard]] auto v30_payload(std::string_view line) -> std::string_view {
    if (!line.starts_with(v30_prefix)) {
        throw std::runtime_error{"expected V3000 record"};
    }

    return line.substr(v30_prefix.size());
}

[[nodiscard]] auto read_v30_line(std::istream& input, std::size_t& line) -> std::string {
    auto result = common::read_line(input, line, "MOL");

    if (!result.starts_with(v30_prefix)) {
        throw std::runtime_error{"expected V3000 record"};
    }

    while (!result.empty() && result.back() == '-') {
        result.pop_back();
        const auto continuation = common::read_line(input, line, "MOL");

        if (!continuation.starts_with(v30_prefix)) {
            throw std::runtime_error{"invalid V3000 continuation"};
        }

        result += continuation.substr(v30_prefix.size());
    }

    return result;
}

[[nodiscard]] auto v30_tokens(std::string_view line) -> std::vector<std::string> {
    std::istringstream input{std::string{v30_payload(line)}};
    std::vector<std::string> result;
    std::string token;

    while (input >> token) {
        result.push_back(std::move(token));
    }

    return result;
}

auto parse_v3000(std::istream& input, std::size_t& line) -> ParsedMolecule {
    if (v30_payload(read_v30_line(input, line)) != "BEGIN CTAB") {
        throw std::runtime_error{"expected V3000 BEGIN CTAB"};
    }

    const auto counts = v30_tokens(read_v30_line(input, line));

    if (counts.size() < 3 || counts[0] != "COUNTS") {
        throw std::runtime_error{"expected V3000 COUNTS"};
    }

    const auto atom_count = common::parse_int(counts[1], "V3000 atom count");
    const auto bond_count = common::parse_int(counts[2], "V3000 bond count");

    if (atom_count <= 0 || bond_count < 0) {
        throw std::runtime_error{"invalid V3000 counts"};
    }

    if (v30_payload(read_v30_line(input, line)) != "BEGIN ATOM") {
        throw std::runtime_error{"expected V3000 BEGIN ATOM"};
    }

    auto result = ParsedMolecule{};
    result.source_connectivity =
        bond_count == 0 ? SourceConnectivity::explicitly_empty : SourceConnectivity::present;
    std::unordered_map<int, std::size_t> atom_indices;
    std::unordered_set<int> bond_ids;
    std::unordered_set<int> ignored_bond_orders;

    for (int index = 0; index < atom_count; ++index) {
        const auto atom = v30_tokens(read_v30_line(input, line));

        if (atom.size() < 6) {
            throw std::runtime_error{"invalid V3000 atom record"};
        }

        const auto source_id = common::parse_int(atom[0], "V3000 atom id");

        if (source_id <= 0 || atom_indices.contains(source_id)) {
            throw std::runtime_error{"duplicate or invalid V3000 atom id"};
        }

        int formal_charge = 0;

        for (std::size_t attribute_index = 6; attribute_index < atom.size(); ++attribute_index) {
            const auto& attribute = atom[attribute_index];

            if (attribute.starts_with("CHG=")) {
                formal_charge =
                    common::parse_int(std::string_view{attribute}.substr(4), "V3000 CHG value");
            }
        }

        atom_indices.emplace(source_id, result.atoms.size());
        result.atoms.push_back(atom_from_symbol(atom[1], formal_charge));
        result.atom_references.push_back(
            SourceAtomReference{.position = static_cast<std::size_t>(index),
                                .id = atom[0],
                                .structural_labels = nullptr});
        result.positions.push_back(
            core::Position{.x = common::parse_double(atom[2], "V3000 x coordinate"),
                           .y = common::parse_double(atom[3], "V3000 y coordinate"),
                           .z = common::parse_double(atom[4], "V3000 z coordinate")});
    }

    if (v30_payload(read_v30_line(input, line)) != "END ATOM") {
        throw std::runtime_error{"expected V3000 END ATOM"};
    }

    auto next = read_v30_line(input, line);
    if (v30_payload(next) == "BEGIN BOND") {
        for (int index = 0; index < bond_count; ++index) {
            const auto bond_line = line + 1;
            const auto bond = v30_tokens(read_v30_line(input, line));

            if (bond.size() < 4) {
                throw std::runtime_error{"invalid V3000 bond record"};
            }

            const auto source_id = common::parse_int(bond[0], "V3000 bond id");
            if (source_id <= 0 || !bond_ids.emplace(source_id).second) {
                throw std::runtime_error{"duplicate or invalid V3000 bond id"};
            }

            const auto first_source_id = common::parse_int(bond[2], "V3000 bond atom");
            const auto second_source_id = common::parse_int(bond[3], "V3000 bond atom");
            const auto first = atom_indices.find(first_source_id);
            const auto second = atom_indices.find(second_source_id);

            if (first == atom_indices.end() || second == atom_indices.end()) {
                throw std::runtime_error{"V3000 bond references an unknown atom"};
            }

            const auto order = common::parse_int(bond[1], "V3000 bond order");
            if (order == 9 || order == 10) {
                if (ignored_bond_orders.emplace(order).second) {
                    result.issues.push_back(ImportIssue{
                        .kind = ImportIssueKind::v3000_bond_order_ignored,
                        .message =
                            "V3000 bond order " + std::to_string(order) +
                            " was ignored because it is not represented in charge calculations; "
                            "the imported graph omits this bond.",
                        .line = bond_line});
                }
                continue;
            }

            result.bonds.emplace_back(first->second, second->second,
                                      core::bond_order_from_value(order == 4 ? 1 : order));
            if (order == 4) {
                result.aromatic_bonds.push_back(result.bonds.size() - 1);
            }
        }

        if (v30_payload(read_v30_line(input, line)) != "END BOND") {
            throw std::runtime_error{"expected V3000 END BOND"};
        }
        next = read_v30_line(input, line);
    } else if (bond_count != 0) {
        throw std::runtime_error{"expected V3000 bond block"};
    }

    // Only the atoms and bonds are read. Other blocks, such as SGROUP and COLLECTION, are ignored,
    // as is anything between the end of the connection table and M  END. A further ATOM or BOND
    // block would change the molecule and is an error.
    while (v30_payload(next) != "END CTAB") {
        if (v30_payload(next) == "BEGIN ATOM" || v30_payload(next) == "BEGIN BOND") {
            throw std::runtime_error{"repeated or misplaced V3000 " +
                                     std::string{v30_payload(next).substr(6)} + " block"};
        }
        next = read_v30_line(input, line);
    }
    while (true) {
        const auto text = common::read_line(input, line, "MOL");
        if (text == "M  END") {
            return result;
        }
        if (text == "$$$$") {
            throw std::runtime_error{"missing M  END"};
        }
    }
}

// The header's program line declares a 2D drawing with "2D" in columns 21-22, for example
// "  -OEChem-10272610592D". A flat 3D molecule with z = 0 does not declare it.
[[nodiscard]] auto declares_2d(const std::string_view program_line) -> bool {
    constexpr std::size_t offset = 20;
    return program_line.size() >= offset + 2 && program_line.substr(offset, 2) == "2D";
}

// A 2D drawing has no geometry, so its record is imported without coordinates; methods that need
// geometry then refuse it, and the others still apply.
[[nodiscard]] auto make_record(ParsedMolecule parsed, MoleculeRecordIdentity identity,
                               const bool drawing, const std::size_t record_start)
    -> ImportedMoleculeRecord {
    auto conformers = std::vector<core::Conformer>{};
    auto conformer_references = std::vector<SourceConformerReference>{};
    if (!parsed.aromatic_bonds.empty()) {
        parsed.issues.push_back(adapters::aromatic_bonds_issue(std::move(parsed.aromatic_bonds)));
    }
    if (drawing) {
        parsed.issues.push_back(ImportIssue{
            .kind = ImportIssueKind::coordinates_2d_ignored,
            .message = "the record is a 2D drawing; its coordinates were not used, so methods "
                       "that need 3D geometry do not apply",
            .line = record_start + 2});
    } else {
        conformers.emplace_back(std::move(parsed.positions), "input");
        conformer_references.push_back(SourceConformerReference{
            .position = 0, .id = std::nullopt, .sites = parsed.atom_references});
    }
    auto metadata = MoleculeImportMetadata{
        .format = MolecularSourceFormat::mol,
        .atoms = std::move(parsed.atom_references),
        .conformers = std::move(conformer_references),
        .policy = {.selection = std::nullopt,
                   .alternate_locations = std::nullopt,
                   .conformers = "all",
                   .bonds = std::nullopt},
        .source_connectivity = parsed.source_connectivity,
        .components = {},
    };
    return adapters::make_record(std::move(parsed.atoms), std::move(parsed.bonds),
                                 std::move(conformers), std::move(identity), {},
                                 std::move(parsed.issues), std::move(metadata));
}

} // namespace

auto parse_mol(std::istream& input, MoleculeRecordIdentity identity) -> ImportedMoleculeRecord {
    std::size_t line = 0;
    return parse_mol(input, std::move(identity), line);
}

auto parse_mol(std::istream& input, MoleculeRecordIdentity identity, std::size_t& line)
    -> ImportedMoleculeRecord {
    const auto record_start = line;
    try {
        const auto name = common::read_line(input, line, "MOL");
        if (identity.record_id.empty()) {
            identity.record_id = name;
        }
        const auto drawing = declares_2d(common::read_line(input, line, "MOL"));
        static_cast<void>(common::read_line(input, line, "MOL"));
        const auto counts = common::read_line(input, line, "MOL");

        ParsedMolecule parsed;

        if (counts.contains("V2000")) {
            parsed = parse_v2000(input, counts, line);
        } else if (counts.contains("V3000")) {
            parsed = parse_v3000(input, line);
        } else {
            throw std::runtime_error{"unsupported MOL version"};
        }

        return make_record(std::move(parsed), identity, drawing, record_start);
    } catch (...) {
        common::rethrow_with_context(identity, line);
    }
}

MolReader::MolReader(std::istream& input, std::string source_name)
    : input_{std::addressof(input)}, source_{std::move(source_name)} {}

auto MolReader::next() -> std::optional<ImportedMoleculeRecord> {
    if (consumed_ || input_->peek() == std::char_traits<char>::eof()) {
        return std::nullopt;
    }

    consumed_ = true;
    auto line = std::size_t{};
    auto result =
        parse_mol(*input_, {.source_name = source_, .record_index = 0, .record_id = {}}, line);
    // Content after M  END is not read as a molecule but must still be ASCII text.
    try {
        while (common::try_read_line(*input_, line).has_value()) {
        }
    } catch (...) {
        common::rethrow_with_context(result.identity, line);
    }
    return result;
}

} // namespace chargefw::adapters::native::mol_input
