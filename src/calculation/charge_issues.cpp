#include <chargefw/calculation/charge_issues.h>

#include "calculation/charge_checks.h"
#include "core/diagnostic_description.h"

#include <cmath>
#include <format>
#include <limits>
#include <stdexcept>
#include <string>

namespace chargefw::calculation {

auto to_string(const ChargeIssueKind value) -> std::string_view {
    switch (value) {
    case ChargeIssueKind::charge_out_of_range:
        return "charge_out_of_range";
    }
    throw std::invalid_argument{"unknown charge issue kind"};
}

namespace detail {

auto find_charge_issues(const core::MoleculeCollection& molecules,
                        const charges::ChargeSet& charges,
                        const FixedChargePartition* fixed_charge_partition)
    -> std::vector<ChargeIssue> {
    auto result = std::vector<ChargeIssue>{};
    for (const auto& assignment : charges.assignments()) {
        const auto molecule_index = assignment.target.molecule_index;
        const auto& molecule = molecules[molecule_index];
        const auto values = assignment.charges.values();
        auto fixed = std::vector<bool>(molecule.atom_count(), false);
        if (fixed_charge_partition != nullptr) {
            for (const auto& source : fixed_charge_partition->targets[molecule_index].sources) {
                fixed[source.atom_index] = true;
            }
        }

        auto atoms = std::vector<std::size_t>{};
        auto furthest = std::size_t{0};
        auto furthest_deviation = 0.0;
        for (std::size_t atom = 0; atom < values.size(); ++atom) {
            if (fixed[atom]) {
                continue;
            }
            // A non-finite charge is as far as possible.
            auto deviation = std::abs(values[atom] - molecule.atom(atom).formal_charge());
            if (std::isnan(deviation)) {
                deviation = std::numeric_limits<double>::infinity();
            }
            if (deviation <= charge_deviation_limit) {
                continue;
            }
            if (atoms.empty() || deviation > furthest_deviation) {
                furthest = atom;
                furthest_deviation = deviation;
            }
            atoms.push_back(atom);
        }
        if (atoms.empty()) {
            continue;
        }

        const auto conformer = assignment.target.conformer_index;
        const auto atom =
            std::format("{} at {:+.2f} e", core::detail::atom_description(molecule, furthest),
                        values[furthest]);
        const auto finding =
            atoms.size() == 1
                ? std::format("{} is more than {} e away from its formal charge", atom,
                              charge_deviation_limit)
                : std::format("{} atoms have charges more than {} e away from their formal "
                              "charges; the furthest is {}",
                              atoms.size(), charge_deviation_limit, atom);
        result.push_back(ChargeIssue{
            .kind = ChargeIssueKind::charge_out_of_range,
            .message =
                std::format("{}{}; the charges are probably wrong, so check the input or "
                            "use another method or parameter set",
                            conformer.has_value() ? std::format("conformer {}: ", *conformer + 1)
                                                  : std::string{},
                            finding),
            .molecule_index = molecule_index,
            .location = {.conformer_index = conformer,
                         .atom_index = furthest,
                         .atom_indices = std::move(atoms)}});
    }
    return result;
}

} // namespace detail
} // namespace chargefw::calculation
