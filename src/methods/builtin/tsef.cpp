#include "methods/builtin/tsef.h"
#include "methods/builtin/dense_solve.h"

#include "features/topology_helpers.h"

#include <chargefw/core/molecule.h>
#include <chargefw/parameters/models/parameter_view.h>

#include <Eigen/Core>

#include <array>
#include <cstddef>
#include <utility>
#include <vector>

namespace chargefw::methods::builtin {
namespace {

constexpr auto alpha = 14.4;

[[nodiscard]] auto k(const int distance) -> double {
    constexpr std::array values{0.556, 0.778, 1.000, 1.053, 1.087, 1.091};

    if (std::cmp_greater(distance, values.size())) {
        return values.back();
    }

    return values.at(static_cast<std::size_t>(distance - 1));
}

[[nodiscard]] auto interaction(const int bond_distance) -> double {
    return alpha * k(bond_distance) / (0.84 * static_cast<double>(bond_distance) + 0.46);
}

} // namespace

auto TSEFMethod::add_method_specific_prerequisite_issues(const MethodPrerequisiteInput& input,
                                                         PrerequisiteResult& result) const -> void {
    if (!features::is_connected(input.prepared_molecule.topology().adjacency())) {
        result.issues.push_back(PrerequisiteIssue{
            .kind = PrerequisiteIssueKind::unsupported_molecule,
            .message = "TSEF requires all atoms to be joined by bonds into one bonded part"});
    }
}

auto TSEFMethod::calculate(const CalculationInput& input) const -> charges::AtomicCharges {
    const auto& molecule = input.molecule();
    const auto& parameters = input.parameters();

    const auto atom_count = molecule.atom_count();

    if (atom_count == 0) {
        return charges::AtomicCharges{std::vector<double>{}};
    }

    const auto electronegativity = parameters.atom("electronegativity");
    const auto hardness = parameters.atom("hardness");

    const auto n = static_cast<Eigen::Index>(atom_count);

    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n + 1, n + 1);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(n + 1);

    const auto bond_distances = features::all_pairs_bond_distances(input.topology().adjacency());

    for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
        const auto i = static_cast<Eigen::Index>(atom_index);

        A(i, i) = hardness[atom_index];
        b(i) = -electronegativity[atom_index];

        for (std::size_t other_atom_index = atom_index + 1; other_atom_index < atom_count;
             ++other_atom_index) {
            const auto bond_distance = bond_distances[atom_index][other_atom_index];
            const auto j = static_cast<Eigen::Index>(other_atom_index);
            const auto value = interaction(bond_distance);

            A(i, j) = value;
            A(j, i) = value;
        }
    }

    const Eigen::VectorXd q = solve_charge_constrained(A, b, input.target_charge()).head(n);

    return charges::AtomicCharges{std::vector<double>{q.data(), q.data() + q.size()}};
}

auto TSEFMethod::estimate(const TargetSize& target, const MethodOptions& /*options*/) const
    -> ResourceEstimate {
    const auto n = target.atoms;
    const auto size = n + 1.0;
    // Dense (n + 1) x (n + 1) system plus all-pairs integer bond distances from one BFS per atom.
    return {.memory_bytes = 8.0 * size * size + 4.0 * n * n,
            .operations = n * (n + target.bonds) + 2.0 / 3.0 * size * size * size};
}

} // namespace chargefw::methods::builtin
