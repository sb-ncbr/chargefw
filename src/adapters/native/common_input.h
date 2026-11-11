#pragma once

#include <chargefw/adapters/molecule_record.h>
#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/conformer.h>

#include "bond_format.h"

#include <cstddef>
#include <istream>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace chargefw::adapters::native::common_input {

// Called from a parser's catch block; preserves standard input exception categories.
[[noreturn]] auto rethrow_with_context(const MoleculeRecordIdentity& identity, std::size_t line)
    -> void;

[[nodiscard]] auto try_read_line(std::istream& input, std::size_t& line)
    -> std::optional<std::string>;

[[nodiscard]] auto read_line(std::istream& input, std::size_t& line, std::string_view record_name)
    -> std::string;

[[nodiscard]] auto trim(std::string_view value) -> std::string_view;

[[nodiscard]] auto parse_int(std::string_view value, std::string_view field) -> int;

[[nodiscard]] auto parse_trimmed_int(std::string_view value, std::string_view field) -> int;

[[nodiscard]] auto parse_double(std::string_view value, std::string_view field) -> double;

[[nodiscard]] auto fixed_field(std::string_view line, std::size_t offset, std::size_t width,
                               std::string_view field) -> std::string_view;

[[nodiscard]] auto bond_order(std::string_view value,
                              ::chargefw::adapters::native::BondFormat format) -> core::BondOrder;

// The warning that a record's aromatic bonds, given by their indices in the record in bond order,
// were read as single bonds.
[[nodiscard]] auto aromatic_bonds_issue(std::vector<std::size_t> bonds) -> ImportIssue;

[[nodiscard]] auto make_record(std::vector<core::Atom> atoms, std::vector<core::Bond> bonds,
                               std::vector<core::Conformer> conformers,
                               MoleculeRecordIdentity identity, std::string name = {},
                               std::vector<ImportIssue> issues = {},
                               std::optional<MoleculeImportMetadata> import_metadata = std::nullopt)
    -> ImportedMoleculeRecord;

} // namespace chargefw::adapters::native::common_input
