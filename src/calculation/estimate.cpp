#include "calculation/assessment_state.h"
#include "calculation/cover_execution.h"
#include "calculation/parallel_for.h"
#include "calculation/target_execution.h"

#include <chargefw/calculation/calculation.h>
#include <chargefw/methods/method.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <numbers>
#include <numeric>
#include <vector>

namespace chargefw::calculation {
namespace {

// Packed biomolecular structures, hydrogens included, hold about 0.1 atoms per cubic angstrom.
constexpr double atom_density = 0.1;
// Output charges retained for every atom of every target.
constexpr double charge_bytes_per_atom = 8.0;

[[nodiscard]] auto atoms_within(const double radius, const double atom_count) -> double {
    const auto sphere_atoms = atom_density * 4.0 / 3.0 * std::numbers::pi * std::pow(radius, 3);
    return std::min(atom_count, std::max(1.0, std::round(sphere_atoms)));
}

} // namespace

auto estimate(const ExecutionPlan& plan, const std::size_t max_threads)
    -> methods::ResourceEstimate {
    const auto& selected = *plan.candidate_;
    const auto& policy = plan.policy();
    const auto reduced = policy.mode() != ExecutionMode::full;
    const auto& molecules = plan.state_->prepared_molecules;
    const auto targets = detail::make_calculation_targets(
        molecules, reduced || selected.method->requirements().coordinates);
    const auto threads = static_cast<double>(detail::thread_count(max_threads));

    auto result = methods::ResourceEstimate{};
    auto full_workspaces = std::vector<double>{};
    auto workspace_memory = 0.0;

    for (const auto& target : targets) {
        const auto& molecule = molecules[target.molecule_index].molecule();
        const auto atoms = static_cast<double>(molecule.atom_count());
        const auto bonds = static_cast<double>(molecule.bond_count());
        result.memory_bytes += charge_bytes_per_atom * atoms;

        if (!reduced) {
            const auto cost = selected.method->estimate({.atoms = atoms, .bonds = bonds},
                                                        selected.method_options);
            result.operations += cost.operations;
            full_workspaces.push_back(cost.memory_bytes);
            continue;
        }

        // Reduced targets run serially, each with one fragment per atom (cutoff) or per cover
        // pivot in parallel; fragments follow the atom density and the target's bonds per atom.
        const auto fragment_atoms = atoms_within(*policy.radius(), atoms);
        const auto fragment_bonds = atoms > 0.0 ? bonds * fragment_atoms / atoms : 0.0;
        const auto fragments =
            policy.mode() == ExecutionMode::cutoff
                ? atoms
                : std::ceil(atoms /
                            std::max(1.0, atoms_within(detail::cover_retained_radius, atoms)));
        const auto cost = selected.method->estimate(
            {.atoms = fragment_atoms, .bonds = fragment_bonds}, selected.method_options);
        result.operations += fragments * cost.operations;
        workspace_memory =
            std::max(workspace_memory, std::min(threads, fragments) * cost.memory_bytes);
    }

    if (!reduced) {
        // Whole targets run in parallel, so the largest workspaces may coincide.
        const auto concurrent = std::min(static_cast<std::size_t>(threads), full_workspaces.size());
        const auto concurrent_end =
            full_workspaces.begin() + static_cast<std::ptrdiff_t>(concurrent);
        std::ranges::partial_sort(full_workspaces, concurrent_end, std::greater{});
        workspace_memory = std::accumulate(full_workspaces.begin(), concurrent_end, 0.0);
    }

    result.memory_bytes += workspace_memory;
    return result;
}

} // namespace chargefw::calculation
