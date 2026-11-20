#pragma once

#include <chargefw/io/molecule_record.h>

#include <charconv>
#include <cstddef>
#include <limits>
#include <optional>
#include <string>
#include <system_error>

#include <gemmi/cifdoc.hpp>

namespace chargefw::io::gemmi::mmcif_labels {

constexpr auto model_id_column = std::size_t{12};

[[nodiscard]] inline auto source_atom_sites(::gemmi::cif::Block& block) -> ::gemmi::cif::Table {
    return block.find("_atom_site.", {"id", "?label_atom_id", "?label_comp_id", "?label_asym_id",
                                      "?label_seq_id", "?auth_atom_id", "?auth_comp_id",
                                      "?auth_asym_id", "?auth_seq_id", "?pdbx_PDB_ins_code",
                                      "?label_alt_id", "?label_entity_id", "?pdbx_PDB_model_num"});
}

[[nodiscard]] inline auto source_value(const ::gemmi::cif::Table::Row& row, const std::size_t index)
    -> std::optional<std::string> {
    if (!row.has(index)) {
        return std::nullopt;
    }
    const auto& token = row[index];
    if (token == "." || token == "?") {
        return token;
    }
    return ::gemmi::cif::as_string(token);
}

[[nodiscard]] inline auto decode_source_labels(const ::gemmi::cif::Table::Row& row)
    -> SourceStructuralLabels {
    return SourceStructuralLabels{.author = SourceHierarchyLabels{.atom = source_value(row, 5),
                                                                  .residue = source_value(row, 6),
                                                                  .chain = source_value(row, 7),
                                                                  .sequence = source_value(row, 8)},
                                  .label = SourceHierarchyLabels{.atom = source_value(row, 1),
                                                                 .residue = source_value(row, 2),
                                                                 .chain = source_value(row, 3),
                                                                 .sequence = source_value(row, 4)},
                                  .entity = source_value(row, 11),
                                  .insertion_code = source_value(row, 9),
                                  .alternate_location = source_value(row, 10),
                                  .segment = std::nullopt};
}

// The formal charge of an _atom_site.pdbx_formal_charge value: zero for a missing or unknown one,
// and nothing for text that is not an integer from -128 to 127, the charges Gemmi can hold. The
// conversion checks for overflow, which Gemmi's own does not.
[[nodiscard]] inline auto formal_charge(const std::string& value) -> std::optional<int> {
    if (::gemmi::cif::is_null(value)) {
        return 0;
    }
    const auto text = ::gemmi::cif::as_string(value);
    const auto* first = text.data() + (text.starts_with('+') ? 1 : 0);
    const auto* last = text.data() + text.size();
    auto charge = 0;
    const auto [end, error] = std::from_chars(first, last, charge);
    // from_chars() reads a minus sign itself, so "+-1" would pass as -1.
    if (error != std::errc{} || end != last || first == last || text.starts_with("+-") ||
        charge < std::numeric_limits<signed char>::min() ||
        charge > std::numeric_limits<signed char>::max()) {
        return std::nullopt;
    }
    return charge;
}

} // namespace chargefw::io::gemmi::mmcif_labels
