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

[[nodiscard]] auto element_symbol(const std::string_view atom_type) -> std::string_view {
    const auto separator = atom_type.find('.');
    const auto symbol = atom_type.substr(0, separator);

    if (symbol.empty() || symbol == "Du" || symbol == "LP" || symbol == "Any" || symbol == "Hal" ||
        symbol == "Het" || symbol == "Hev") {
        throw std::runtime_error{"unsupported MOL2 atom type '" + std::string{atom_type} + "'"};
    }

    return symbol;
}

auto read_to_atom_section(std::istream& input, std::size_t& line) -> void {
    while (true) {
        const auto section = common::read_line(input, line, "MOL2");
        if (section == atom_marker) {
            return;
        }
        if (section == molecule_marker) {
            throw std::runtime_error{"missing MOL2 ATOM section"};
        }
        // Skip unrelated MOL2 sections: this adapter reads only MOLECULE, ATOM, and BOND data.
    }
}

[[nodiscard]] auto parse_record(std::istream& input, MoleculeRecordIdentity identity)
    -> ImportedMoleculeRecord {
    std::size_t line = 1;

    const auto name = common::read_line(input, line, "MOL2");
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

    read_to_atom_section(input, line);

    std::vector<core::Atom> atoms;
    std::vector<core::Position> positions;
    std::vector<SourceAtomReference> atom_references;
    std::unordered_map<int, std::size_t> atom_indices;
    bool has_partial_charges = false;
    atoms.reserve(static_cast<std::size_t>(atom_count));
    positions.reserve(static_cast<std::size_t>(atom_count));
    atom_references.reserve(static_cast<std::size_t>(atom_count));

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

        if (source_id <= 0 || atom_indices.contains(source_id)) {
            throw std::runtime_error{"duplicate or invalid MOL2 atom id"};
        }

        if (!charge_text.empty()) {
            const auto charge = common::parse_double(charge_text, "MOL2 atom charge");

            if (charge != 0.0) {
                has_partial_charges = true;
            }
        }

        const auto& element = core::periodic_table().element(element_symbol(atom_type));
        atom_indices.emplace(source_id, atoms.size());
        atoms.emplace_back(element.atomic_number, 0, std::move(atom_name));
        atom_references.push_back(SourceAtomReference{.position = static_cast<std::size_t>(index),
                                                      .id = std::move(source_id_text),
                                                      .structural_labels = nullptr});
        positions.push_back(core::Position{.x = common::parse_double(x_text, "MOL2 x coordinate"),
                                           .y = common::parse_double(y_text, "MOL2 y coordinate"),
                                           .z = common::parse_double(z_text, "MOL2 z coordinate")});
    }

    if (common::read_line(input, line, "MOL2") != bond_marker) {
        throw std::runtime_error{"expected MOL2 BOND section"};
    }

    std::vector<core::Bond> bonds;
    bonds.reserve(static_cast<std::size_t>(bond_count));

    for (int index = 0; index < bond_count; ++index) {
        std::istringstream bond_line{common::read_line(input, line, "MOL2")};
        std::string bond_id;
        std::string first_text;
        std::string second_text;
        std::string type;

        if (!(bond_line >> bond_id >> first_text >> second_text >> type)) {
            throw std::runtime_error{"invalid MOL2 bond record"};
        }

        const auto first_source_id = common::parse_int(first_text, "MOL2 bond atom");
        const auto second_source_id = common::parse_int(second_text, "MOL2 bond atom");
        const auto first = atom_indices.find(first_source_id);
        const auto second = atom_indices.find(second_source_id);

        if (first == atom_indices.end() || second == atom_indices.end()) {
            throw std::runtime_error{"MOL2 bond references an unknown atom"};
        }

        bonds.emplace_back(first->second, second->second,
                           common::bond_order(type, BondFormat::mol2));
    }

    if (identity.record_id.empty()) {
        identity.record_id = name;
    }

    auto diagnostics = std::vector<MoleculeRecordDiagnostic>{};

    if (has_partial_charges) {
        diagnostics.push_back(MoleculeRecordDiagnostic{
            .code = "partial_charges_ignored",
            .message = "MOL2 partial charges were ignored; they are not formal charges and "
                       "were not used for calculation.",
            .line = std::nullopt});
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
        .source_connectivity =
            bond_count == 0 ? SourceConnectivity::explicitly_empty : SourceConnectivity::present,
        .components = {},
    };
    return common::make_record(
        std::move(atoms), std::move(bonds), {core::Conformer{std::move(positions), "input"}},
        std::move(identity), {}, std::move(diagnostics), std::move(metadata));
}

} // namespace

Mol2Reader::Mol2Reader(std::istream& input, std::string source)
    : input_{std::addressof(input)}, source_{std::move(source)} {}

auto Mol2Reader::next() -> std::optional<ImportedMoleculeRecord> {
    std::string line;
    while (std::getline(*input_, line)) {
        if (line.ends_with('\r')) {
            line.pop_back();
        }
        if (line == molecule_marker) {
            return parse_record(*input_, MoleculeRecordIdentity{.source = source_,
                                                                .record_index = record_index_++,
                                                                .record_id = {}});
        }
    }

    return std::nullopt;
}

} // namespace chargefw::adapters::native::mol2_input
