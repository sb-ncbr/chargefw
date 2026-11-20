#pragma once

#include <chargefw/io/molecule_record.h>
#include <chargefw/io/parse_error.h>

#include <gemmi/model.hpp>

#include <charconv>
#include <cstdlib>
#include <optional>
#include <string>
#include <string_view>
#include <system_error>
#include <tuple>

// How the atoms of a PDB file are identified by their labels, for the reader and for the
// conversion of the file to mmCIF, which must find the same atoms again.
namespace chargefw::io::gemmi::pdb_labels {

[[nodiscard]] inline auto optional_text(const std::string_view value)
    -> std::optional<std::string> {
    return value.empty() ? std::nullopt : std::optional<std::string>{value};
}

[[nodiscard]] inline auto optional_character(const char value) -> std::optional<std::string> {
    return value == '\0' || value == ' ' ? std::nullopt
                                         : std::optional<std::string>{std::string{value}};
}

[[nodiscard]] inline auto normalized_sequence(const std::optional<std::string>& sequence)
    -> std::optional<std::string> {
    if (!sequence.has_value()) {
        return std::nullopt;
    }
    if (sequence->front() < 'A') {
        auto value = int{};
        const auto [end, error] =
            std::from_chars(sequence->data(), sequence->data() + sequence->size(), value);
        if (error != std::errc{} || end != sequence->data() + sequence->size()) {
            throw ParseError{ParseErrorKind::malformed,
                             "PDB source contains an invalid residue sequence ID",
                             {}};
        }
        return std::to_string(value);
    }

    const auto value = std::strtol(sequence->c_str(), nullptr, 36) - 466560L + 10000L;
    return std::to_string(value);
}

struct StructuralLabelsLess {
    [[nodiscard]] auto operator()(const SourceStructuralLabels& lhs,
                                  const SourceStructuralLabels& rhs) const -> bool {
        return std::tie(lhs.author.atom, lhs.author.residue, lhs.author.chain, lhs.author.sequence,
                        lhs.label.atom, lhs.label.residue, lhs.label.chain, lhs.label.sequence,
                        lhs.entity, lhs.insertion_code, lhs.alternate_location, lhs.segment) <
               std::tie(rhs.author.atom, rhs.author.residue, rhs.author.chain, rhs.author.sequence,
                        rhs.label.atom, rhs.label.residue, rhs.label.chain, rhs.label.sequence,
                        rhs.entity, rhs.insertion_code, rhs.alternate_location, rhs.segment);
    }
};

[[nodiscard]] inline auto selected_labels(const ::gemmi::const_CRA& site)
    -> SourceStructuralLabels {
    return SourceStructuralLabels{
        .author = SourceHierarchyLabels{.atom = optional_text(site.atom->name),
                                        .residue = optional_text(site.residue->name),
                                        .chain = optional_text(site.chain->name),
                                        .sequence = optional_text(site.residue->seqid.num.str())},
        .label = {},
        .entity = std::nullopt,
        .insertion_code = optional_character(site.residue->seqid.icode),
        .alternate_location = optional_character(site.atom->altloc),
        .segment = optional_text(site.residue->segment)};
}

} // namespace chargefw::io::gemmi::pdb_labels
