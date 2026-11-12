#pragma once

#include <chargefw/core/issue_location.h>
#include <chargefw/core/molecule.h>

#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace chargefw::adapters {

class PortableId {
  public:
    using Value = std::variant<std::int64_t, std::string>;

    PortableId() = default;
    // Implicit on purpose: an ID is written as the integer or string it holds, and every integer or
    // string is a valid ID.
    // NOLINTBEGIN(cppcoreguidelines-explicit-constructor)
    PortableId(std::int64_t value) : value_{value} {}
    PortableId(std::string value) : value_{std::move(value)} {}
    PortableId(const char* value) : value_{std::string{value}} {}
    // NOLINTEND(cppcoreguidelines-explicit-constructor)

    [[nodiscard]] auto empty() const noexcept -> bool {
        return !value_.has_value();
    }
    [[nodiscard]] auto value() const noexcept -> const std::optional<Value>& {
        return value_;
    }
    [[nodiscard]] auto display_string() const -> std::string {
        if (!value_.has_value()) {
            return {};
        }
        return std::visit(
            [](const auto& value) -> std::string {
                using T = std::decay_t<decltype(value)>;
                if constexpr (std::is_same_v<T, std::string>) {
                    return value;
                } else {
                    return std::to_string(value);
                }
            },
            *value_);
    }

    [[nodiscard]] auto operator==(const PortableId&) const -> bool = default;
    [[nodiscard]] friend auto operator==(const PortableId& id, const std::string_view value)
        -> bool {
        return id.value_.has_value() && std::holds_alternative<std::string>(*id.value_) &&
               std::get<std::string>(*id.value_) == value;
    }
    [[nodiscard]] friend auto operator==(const std::string_view value, const PortableId& id)
        -> bool {
        return id == value;
    }
    [[nodiscard]] friend auto operator==(const PortableId& id, const char* value) -> bool {
        return id == std::string_view{value};
    }
    [[nodiscard]] friend auto operator==(const char* value, const PortableId& id) -> bool {
        return id == std::string_view{value};
    }
    [[nodiscard]] friend auto operator==(const PortableId& id, const std::int64_t value) -> bool {
        return id.value_.has_value() && std::holds_alternative<std::int64_t>(*id.value_) &&
               std::get<std::int64_t>(*id.value_) == value;
    }

  private:
    std::optional<Value> value_;
};

enum class MolecularSourceFormat : std::uint8_t {
    molecule_json,
    mol,
    sdf,
    mol2,
    pdb,
    mmcif,
};

// The name of the format in result JSON and in Python, such as "molecule-json".
[[nodiscard]] auto to_string(MolecularSourceFormat value) -> std::string_view;

// The format with the given name. Throws std::invalid_argument for an unknown name.
[[nodiscard]] auto molecular_source_format_from_string(std::string_view value)
    -> MolecularSourceFormat;

enum class SourceConnectivity : std::uint8_t {
    absent,
    explicitly_empty,
    present,
};

// The name of the connectivity state in result JSON and in Python, such as "explicitly-empty".
[[nodiscard]] auto to_string(SourceConnectivity value) -> std::string_view;

struct SourceHierarchyLabels {
    std::optional<std::string> atom;
    std::optional<std::string> residue;
    std::optional<std::string> chain;
    std::optional<std::string> sequence;

    [[nodiscard]] auto operator==(const SourceHierarchyLabels&) const -> bool = default;
};

struct SourceStructuralLabels {
    SourceHierarchyLabels author;
    SourceHierarchyLabels label;
    std::optional<std::string> entity;
    std::optional<std::string> insertion_code;
    std::optional<std::string> alternate_location;
    std::optional<std::string> segment;

    [[nodiscard]] auto operator==(const SourceStructuralLabels&) const -> bool = default;
};

// Source positions are zero-based within one record. IDs retain explicit source tokens when the
// format provides them; position remains authoritative when no portable ID exists.
struct SourceAtomReference {
    std::size_t position = 0;
    std::optional<std::string> id;
    // PDB and mmCIF only. Readers share one immutable object among sites with equal labels, such
    // as an atom's sites in the first and later models; replace it instead of editing it.
    std::shared_ptr<const SourceStructuralLabels> structural_labels;

    // Compares label values, not objects, so records read separately compare equal.
    [[nodiscard]] friend auto operator==(const SourceAtomReference& first,
                                         const SourceAtomReference& second) -> bool {
        const auto same_labels =
            first.structural_labels == second.structural_labels ||
            (first.structural_labels != nullptr && second.structural_labels != nullptr &&
             *first.structural_labels == *second.structural_labels);
        return first.position == second.position && first.id == second.id && same_labels;
    }
};

struct SourceConformerReference {
    std::size_t position = 0;
    std::optional<std::string> id;
    std::vector<SourceAtomReference> sites;

    [[nodiscard]] auto operator==(const SourceConformerReference&) const -> bool = default;
};

// One selected structural residue instance. Its vector position in metadata.components is its
// record-local ID; atom_indices are in calculation-atom order.
struct SourceComponentInstance {
    std::string component_id;
    std::vector<std::size_t> atom_indices;

    [[nodiscard]] auto operator==(const SourceComponentInstance&) const -> bool = default;
};

// Small adapter-owned metadata retained with the normalized molecule. Structural labels are added
// by the structural-reader boundary rather than by core::Molecule.
struct MoleculeImportMetadata {
    MolecularSourceFormat format = MolecularSourceFormat::molecule_json;
    std::vector<SourceAtomReference> atoms;
    std::vector<SourceConformerReference> conformers;
    std::optional<std::string> record_selection;
    std::optional<std::string> alternate_location_selection;
    std::optional<std::string> conformer_selection;
    std::optional<std::string> bond_strategy;
    SourceConnectivity source_connectivity = SourceConnectivity::absent;
    std::vector<SourceComponentInstance> components;

    [[nodiscard]] auto operator==(const MoleculeImportMetadata&) const -> bool = default;
};

// Identifies one source record without imposing a file-format or toolkit dependency. The source
// and a text record ID must be UTF-8; result JSON cannot hold other text.
struct MoleculeRecordIdentity {
    std::string source;
    std::size_t record_index = 0;
    PortableId record_id;

    [[nodiscard]] auto operator==(const MoleculeRecordIdentity&) const -> bool = default;
};

enum class ImportIssueKind : std::uint8_t {
    partial_charges_ignored,
    coordinates_2d_ignored,
    v3000_bond_order_ignored,
    aromatic_bonds_converted,
    symmetry_connections_ignored,
    unbonded_residue_atoms,
};

// The code of the kind in result JSON, such as "partial_charges_ignored".
[[nodiscard]] auto to_string(ImportIssueKind value) -> std::string_view;

// A caution that a reader attached to a record, such as file data it did not use.
struct ImportIssue {
    ImportIssueKind kind;
    std::string message;
    std::optional<std::size_t> line; // One-based line in the source file.
    // The atoms or bonds of the record the issue concerns, such as the aromatic bonds read as
    // single.
    core::IssueLocation location = {};

    [[nodiscard]] auto operator==(const ImportIssue&) const -> bool = default;
};

// A successful import result. Molecule ownership and source identity travel together so output
// adapters can retain source ordering and identity.
struct ImportedMoleculeRecord {
    core::Molecule molecule;
    MoleculeRecordIdentity identity;
    // Explicit caller IDs are separate from source mappings, whose source positions and tokens
    // remain authoritative for imported records.
    std::optional<std::vector<PortableId>> caller_atom_ids;
    std::vector<ImportIssue> issues;
    std::optional<MoleculeImportMetadata> import_metadata;

    [[nodiscard]] auto operator==(const ImportedMoleculeRecord&) const -> bool = default;
};

} // namespace chargefw::adapters
