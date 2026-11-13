#include <chargefw/adapters/native/json_input.h>

#include "common_input.h"

#include "adapters/imported_records.h"

#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/position.h>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <istream>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace chargefw::adapters::native::json_input {
namespace {

namespace common = chargefw::adapters::native::common_input;

using Json = nlohmann::json;

[[nodiscard]] auto type_name(const Json& value) -> std::string {
    return value.type_name();
}

auto require_object(const Json& value, const std::string& context) -> void {
    if (!value.is_object()) {
        throw std::runtime_error{context + " must be an object, got " + type_name(value)};
    }
}

auto require_array(const Json& value, const std::string& context) -> void {
    if (!value.is_array()) {
        throw std::runtime_error{context + " must be an array, got " + type_name(value)};
    }
}

[[nodiscard]] auto member(const Json& object, const char* name, const std::string& context)
    -> const Json& {
    const auto iterator = object.find(name);

    if (iterator == object.end()) {
        throw std::runtime_error{context + " is missing required member '" + name + "'"};
    }

    return *iterator;
}

[[nodiscard]] auto optional_member(const Json& object, const char* name) -> const Json* {
    const auto iterator = object.find(name);
    return iterator == object.end() ? nullptr : std::addressof(*iterator);
}

[[nodiscard]] auto require_string(const Json& value, const std::string& context) -> std::string {
    if (!value.is_string()) {
        throw std::runtime_error{context + " must be a string, got " + type_name(value)};
    }

    return value.get<std::string>();
}

[[nodiscard]] auto require_integer(const Json& value, const std::string& context) -> int {
    if (!value.is_number_integer()) {
        throw std::runtime_error{context + " must be an integer, got " + type_name(value)};
    }

    if (value.is_number_unsigned()) {
        const auto number = value.get<std::uint64_t>();
        if (number > static_cast<std::uint64_t>(std::numeric_limits<int>::max())) {
            throw std::runtime_error{context + " is outside the supported integer range"};
        }
        return static_cast<int>(number);
    }

    const auto number = value.get<std::int64_t>();
    if (number < static_cast<std::int64_t>(std::numeric_limits<int>::min()) ||
        number > static_cast<std::int64_t>(std::numeric_limits<int>::max())) {
        throw std::runtime_error{context + " is outside the supported integer range"};
    }
    return static_cast<int>(number);
}

[[nodiscard]] auto require_index(const Json& value, const std::string& context) -> std::size_t {
    if (!value.is_number_unsigned()) {
        throw std::runtime_error{context + " must be a non-negative integer, got " +
                                 type_name(value)};
    }

    return value.get<std::size_t>();
}

[[nodiscard]] auto require_number(const Json& value, const std::string& context) -> double {
    if (!value.is_number()) {
        throw std::runtime_error{context + " must be a number, got " + type_name(value)};
    }

    const auto result = value.get<double>();

    if (!std::isfinite(result)) {
        throw std::runtime_error{context + " must be finite"};
    }

    return result;
}

[[nodiscard]] auto index_context(const std::string& context, const std::size_t index)
    -> std::string {
    return context + "[" + std::to_string(index) + "]";
}

[[nodiscard]] auto parse_molecule(const Json& value, MoleculeRecordIdentity identity,
                                  const ConformerSelection conformer_selection)
    -> ImportedMoleculeRecord {
    const std::string context{"molecules[" + std::to_string(identity.record_index) + "]"};
    require_object(value, context);

    const auto* id = optional_member(value, "id");
    if (id != nullptr) {
        identity.record_id = require_string(*id, context + ".id");
    }

    auto name = std::string{};
    const auto* name_value = optional_member(value, "name");
    if (name_value != nullptr) {
        name = require_string(*name_value, context + ".name");
    }

    const auto& atoms_value = member(value, "atoms", context);
    require_array(atoms_value, context + ".atoms");
    if (atoms_value.empty()) {
        throw std::runtime_error{context + ".atoms must not be empty"};
    }

    std::vector<core::Atom> atoms;
    std::vector<SourceAtomReference> atom_references;
    atoms.reserve(atoms_value.size());
    atom_references.reserve(atoms_value.size());
    for (std::size_t index = 0; index < atoms_value.size(); ++index) {
        const auto atom_context = index_context(context + ".atoms", index);
        const auto& atom_value = atoms_value[index];
        require_object(atom_value, atom_context);
        const auto atomic_number = require_integer(
            member(atom_value, "atomic_number", atom_context), atom_context + ".atomic_number");
        const auto formal_charge = require_integer(
            member(atom_value, "formal_charge", atom_context), atom_context + ".formal_charge");
        atoms.emplace_back(atomic_number, formal_charge);
        atom_references.push_back(SourceAtomReference{
            .position = index, .id = std::nullopt, .structural_labels = nullptr});
    }

    std::vector<core::Bond> bonds;
    const auto* bonds_value = optional_member(value, "bonds");
    if (bonds_value != nullptr) {
        require_array(*bonds_value, context + ".bonds");
        bonds.reserve(bonds_value->size());
        for (std::size_t index = 0; index < bonds_value->size(); ++index) {
            const auto bond_context = index_context(context + ".bonds", index);
            const auto& bond_value = (*bonds_value)[index];
            require_object(bond_value, bond_context);
            const auto& endpoints = member(bond_value, "atoms", bond_context);
            require_array(endpoints, bond_context + ".atoms");
            if (endpoints.size() != 2) {
                throw std::runtime_error{bond_context + ".atoms must contain exactly two indices"};
            }
            bonds.emplace_back(
                require_index(endpoints[0], bond_context + ".atoms[0]"),
                require_index(endpoints[1], bond_context + ".atoms[1]"),
                core::bond_order_from_value(require_integer(
                    member(bond_value, "order", bond_context), bond_context + ".order")));
        }
    }

    std::vector<core::Conformer> conformers;
    std::vector<SourceConformerReference> conformer_references;
    const auto* conformers_value = optional_member(value, "conformers");
    if (conformers_value != nullptr) {
        require_array(*conformers_value, context + ".conformers");
        const auto conformer_count = conformer_selection == ConformerSelection::first
                                         ? std::min<std::size_t>(1, conformers_value->size())
                                         : conformers_value->size();
        conformers.reserve(conformer_count);
        conformer_references.reserve(conformer_count);
        for (std::size_t index = 0; index < conformer_count; ++index) {
            const auto conformer_context = index_context(context + ".conformers", index);
            const auto& conformer_value = (*conformers_value)[index];
            require_object(conformer_value, conformer_context);
            auto conformer_name = std::string{};
            const auto* conformer_id = optional_member(conformer_value, "id");
            auto source_conformer_id = std::optional<std::string>{};
            if (conformer_id != nullptr) {
                conformer_name = require_string(*conformer_id, conformer_context + ".id");
                source_conformer_id = conformer_name;
            }
            const auto& coordinates = member(conformer_value, "coordinates", conformer_context);
            require_array(coordinates, conformer_context + ".coordinates");
            std::vector<core::Position> positions;
            positions.reserve(coordinates.size());
            for (std::size_t coordinate_index = 0; coordinate_index < coordinates.size();
                 ++coordinate_index) {
                const auto coordinate_context =
                    index_context(conformer_context + ".coordinates", coordinate_index);
                const auto& coordinate = coordinates[coordinate_index];
                require_array(coordinate, coordinate_context);
                if (coordinate.size() != 3) {
                    throw std::runtime_error{coordinate_context +
                                             " must contain exactly three values"};
                }
                positions.push_back(
                    core::Position{.x = require_number(coordinate[0], coordinate_context + "[0]"),
                                   .y = require_number(coordinate[1], coordinate_context + "[1]"),
                                   .z = require_number(coordinate[2], coordinate_context + "[2]")});
            }
            conformers.emplace_back(std::move(positions), std::move(conformer_name));
            conformer_references.push_back(SourceConformerReference{
                .position = index,
                .id = std::move(source_conformer_id),
                .sites = atom_references,
            });
        }
    }

    auto metadata = MoleculeImportMetadata{
        .format = MolecularSourceFormat::molecule_json,
        .atoms = std::move(atom_references),
        .conformers = std::move(conformer_references),
        .record_selection = std::nullopt,
        .alternate_location_selection = std::nullopt,
        .conformer_selection = std::string{to_string(conformer_selection)},
        .bond_strategy = std::nullopt,
        .source_connectivity = bonds_value == nullptr
                                   ? SourceConnectivity::absent
                                   : (bonds_value->empty() ? SourceConnectivity::explicitly_empty
                                                           : SourceConnectivity::present),
        .components = {},
    };
    return adapters::make_record(std::move(atoms), std::move(bonds), std::move(conformers),
                                 std::move(identity), std::move(name), {}, std::move(metadata));
}

} // namespace

struct JsonReader::Impl {
    Impl(std::istream& input, std::string source, ConformerSelection conformers);
    auto next() -> std::optional<ImportedMoleculeRecord>;

    nlohmann::json molecule_records_;
    std::string source_;
    ConformerSelection conformers_;
    std::size_t record_index_ = 0;
};

JsonReader::JsonReader(std::istream& input, std::string source, const ConformerSelection conformers)
    : impl_{std::make_unique<Impl>(input, std::move(source), conformers)} {}

JsonReader::JsonReader(JsonReader&&) noexcept = default;
auto JsonReader::operator=(JsonReader&&) noexcept -> JsonReader& = default;
JsonReader::~JsonReader() = default;

auto JsonReader::next() -> std::optional<ImportedMoleculeRecord> {
    return impl_ ? impl_->next() : std::nullopt;
}

JsonReader::Impl::Impl(std::istream& input, std::string source, const ConformerSelection conformers)
    : source_{std::move(source)}, conformers_{conformers} {
    auto document = Json::parse(input, nullptr, false);
    if (document.is_discarded()) {
        throw std::runtime_error{"invalid JSON document"};
    }
    require_object(document, "document");

    if (require_string(member(document, "schema_version", "document"), "document.schema_version") !=
        "1.0") {
        throw std::runtime_error{"unsupported schema_version; expected '1.0'"};
    }

    const std::string document_context{"document"};
    require_array(member(document, "molecules", document_context), "document.molecules");
    molecule_records_ = std::move(document["molecules"]);
}

auto JsonReader::Impl::next() -> std::optional<ImportedMoleculeRecord> {
    if (record_index_ == molecule_records_.size()) {
        return std::nullopt;
    }

    const auto current_record_index = record_index_++;
    auto record = parse_molecule(molecule_records_[current_record_index],
                                 MoleculeRecordIdentity{.source = source_,
                                                        .record_index = current_record_index,
                                                        .record_id = {}},
                                 conformers_);
    // The reader only moves forward; release the parsed record's document.
    molecule_records_[current_record_index] = nullptr;
    return record;
}

} // namespace chargefw::adapters::native::json_input
