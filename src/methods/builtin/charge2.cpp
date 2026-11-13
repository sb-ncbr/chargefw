#include "methods/builtin/charge2.h"
#include "methods/builtin/prerequisite_helpers.h"

#include <chargefw/core/periodic_table.h>
#include <chargefw/features/topology_features.h>
#include <chargefw/parameters/models/parameter_view.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace chargefw::methods::builtin {
namespace {

constexpr auto hydrogen_atomic_number = 1;
// Fixed by the canonical Charge2 parameter set.
constexpr auto hydrogen_electronegativity_reference = 7.17;
constexpr auto convergence_tolerance = 1.0e-5;

[[nodiscard]] auto alpha_denominator(const core::Atom& first, const core::Atom& second,
                                     const core::PeriodicTable& table, const double a1,
                                     const double a2, const double a3) -> double {
    if (table.element(first.atomic_number()).period == 2 &&
        table.element(second.atomic_number()).period == 2) {
        return a1;
    }
    if (first.atomic_number() == hydrogen_atomic_number ||
        second.atomic_number() == hydrogen_atomic_number) {
        return a2;
    }
    return a3;
}

struct RemoteEffect {
    std::size_t first_neighbor;
    std::size_t substituent;
};

struct Charge2Shells {
    std::vector<std::size_t> one_bond;
    std::vector<RemoteEffect> two_bonds;
    std::vector<RemoteEffect> three_bonds;
};

[[nodiscard]] auto precompute_charge2_shells(const features::TopologyFeatures& topology)
    -> std::vector<Charge2Shells> {
    const auto atom_count = topology.molecule().atom_count();
    auto shells_by_atom = std::vector<Charge2Shells>(atom_count);

    for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
        auto& shells = shells_by_atom[atom_index];
        for (const auto first_neighbor : topology.neighbor_indices(atom_index)) {
            shells.one_bond.push_back(first_neighbor);

            for (const auto second_neighbor : topology.neighbor_indices(first_neighbor)) {
                if (second_neighbor == atom_index) {
                    continue;
                }
                shells.two_bonds.push_back({first_neighbor, second_neighbor});

                for (const auto third_neighbor : topology.neighbor_indices(second_neighbor)) {
                    if (third_neighbor == atom_index || third_neighbor == first_neighbor) {
                        continue;
                    }
                    shells.three_bonds.push_back({first_neighbor, third_neighbor});
                }
            }
        }
    }

    return shells_by_atom;
}

} // namespace

auto Charge2Method::add_method_specific_prerequisite_issues(const MethodPrerequisiteInput& input,
                                                            PrerequisiteResult& result) const
    -> void {
    detail::add_component_neutrality_prerequisite_issue(input, result, "Charge2");
}

auto Charge2Method::calculate(const CalculationInput& input) const -> charges::AtomicCharges {
    const auto iterations = input.method_options().get<int>("iters");

    const auto& molecule = input.molecule();
    const auto& parameters = input.parameters();
    const auto& table = core::periodic_table();

    const auto chi = parameters.atom("chi");
    const auto p0 = parameters.atom("P0");
    const auto q0 = parameters.atom("q0");

    const auto a1 = parameters.common("a1");
    const auto a2 = parameters.common("a2");
    const auto a3 = parameters.common("a3");
    const auto b = parameters.common("b");
    const auto c = parameters.common("c");
    const auto alpha = parameters.common("alpha");

    const auto atom_count = molecule.atom_count();
    const auto shells_by_atom = precompute_charge2_shells(input.topology());

    auto q = std::vector(atom_count, 0.0);
    auto next_q = std::vector(atom_count, 0.0);

    for (auto iteration = 0; iteration < iterations; ++iteration) {
        std::ranges::fill(next_q, 0.0);

        for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
            const auto& shells = shells_by_atom[atom_index];

            for (const auto first_neighbor : shells.one_bond) {
                if (first_neighbor < atom_index) {
                    continue;
                }
                const auto denominator = alpha_denominator(
                    molecule.atom(atom_index), molecule.atom(first_neighbor), table, a1, a2, a3);
                const auto transfer = (chi[first_neighbor] - chi[atom_index]) / denominator;
                next_q[atom_index] += transfer;
                next_q[first_neighbor] -= transfer;
            }

            const auto polarizability =
                std::max(0.0, p0[atom_index] * (1.0 + alpha * (q0[atom_index] - q[atom_index])));

            for (const auto& [first_neighbor, substituent] : shells.two_bonds) {
                const auto transfer =
                    (chi[substituent] - hydrogen_electronegativity_reference) * polarizability / b;
                next_q[atom_index] += transfer;
                next_q[first_neighbor] -= transfer;
            }

            for (const auto& [first_neighbor, substituent] : shells.three_bonds) {
                const auto transfer = (chi[substituent] - hydrogen_electronegativity_reference) *
                                      polarizability / b / c;
                next_q[atom_index] += transfer;
                next_q[first_neighbor] -= transfer;
            }
        }

        auto maximum_change = 0.0;
        for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
            maximum_change = std::max(maximum_change, std::abs(next_q[atom_index] - q[atom_index]));
        }
        q.swap(next_q);
        if (maximum_change <= convergence_tolerance) {
            break;
        }
    }

    return charges::AtomicCharges{std::move(q)};
}

auto Charge2Method::estimate(const TargetSize& target, const MethodOptions& options) const
    -> ResourceEstimate {
    // An atom with at most four bonds has up to 4 one-bond, 12 two-bond, and 36 three-bond paths.
    // Per atom: two charge vectors, three path-vector headers, and the stored paths.
    constexpr auto bytes_per_atom = 16.0 + 72.0 + 8.0 * 4 + 16.0 * (12 + 36);
    constexpr auto paths_per_atom = 4.0 + 12.0 + 36.0;
    const auto iterations = options.get<int>("iters");
    return {.memory_bytes = bytes_per_atom * target.atoms,
            .operations = iterations * paths_per_atom * target.atoms};
}

} // namespace chargefw::methods::builtin
