#include "methods/builtin/abeem.h"
#include "methods/builtin/dense_solve.h"

#include "methods/builtin/prerequisite_helpers.h"

#include <chargefw/core/molecule.h>
#include <chargefw/core/periodic_table.h>
#include <chargefw/core/position.h>
#include <chargefw/parameters/models/parameter_view.h>

#include <Eigen/Core>

#include <cstddef>
#include <vector>

namespace chargefw::methods::builtin {
namespace {

[[nodiscard]] auto contains_atom(const core::Bond& bond, const std::size_t atom_index) noexcept
    -> bool {
    return bond.first_atom_index() == atom_index || bond.second_atom_index() == atom_index;
}

[[nodiscard]] auto weighted_bond_center(const core::Molecule& molecule,
                                        const features::ConformerFeatures& geometry,
                                        const std::size_t bond_index) -> core::Position {
    const auto& bond = molecule.bond(bond_index);
    const auto& first_atom = molecule.atom(bond.first_atom_index());
    const auto& second_atom = molecule.atom(bond.second_atom_index());

    const auto& table = core::periodic_table();
    const auto& first_element = table.element(first_atom.atomic_number());
    const auto& second_element = table.element(second_atom.atomic_number());

    const auto first_weight = first_element.covalent_radius;
    const auto second_weight = second_element.covalent_radius;
    const auto weight_sum = first_weight + second_weight;

    const auto& [x1, y1, z1] = geometry.position(bond.first_atom_index());
    const auto& [x2, y2, z2] = geometry.position(bond.second_atom_index());

    return core::Position{.x = (second_weight * x1 + first_weight * x2) / weight_sum,
                          .y = (second_weight * y1 + first_weight * y2) / weight_sum,
                          .z = (second_weight * z1 + first_weight * z2) / weight_sum};
}

[[nodiscard]] auto atom_bond_distance(const features::ConformerFeatures& geometry,
                                      const std::size_t atom_index,
                                      const core::Position& bond_center) -> double {
    return core::distance(geometry.position(atom_index), bond_center);
}

[[nodiscard]] auto interaction(const double k, const double distance_value) -> double {
    return k / distance_value;
}

} // namespace

auto ABEEMMethod::add_method_specific_prerequisite_issues(const MethodPrerequisiteInput& input,
                                                          PrerequisiteResult& result) const
    -> void {
    detail::add_element_prerequisite_issues(
        input, result, "ABEEM requires a positive covalent radius",
        [](const core::Element& element) -> bool { return element.covalent_radius > 0.0; });
}

auto ABEEMMethod::calculate(const CalculationInput& input) const -> charges::AtomicCharges {
    const auto& molecule = input.molecule();
    const auto& geometry = input.geometry();
    const auto& parameters = input.parameters();

    const auto atom_count = molecule.atom_count();
    const auto bond_count = molecule.bond_count();

    if (atom_count == 0) {
        return charges::AtomicCharges{std::vector<double>{}};
    }

    const auto k = parameters.common("k");

    if (k <= 0.0) {
        throw std::logic_error{"ABEEM common parameter 'k' must be positive"};
    }

    const auto atom_a = parameters.atom("a");
    const auto atom_b = parameters.atom("b");
    const auto atom_c = parameters.atom("c");

    const auto bond_a = parameters.bond("A");
    const auto bond_b = parameters.bond("B");
    const auto bond_c = parameters.bond("C");
    const auto bond_d = parameters.bond("D");

    const auto n = static_cast<Eigen::Index>(atom_count);
    const auto m = static_cast<Eigen::Index>(bond_count);
    const auto system_size = n + m + 1;

    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(system_size, system_size);
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(system_size);

    auto bond_centers = std::vector<core::Position>{};
    bond_centers.reserve(bond_count);

    for (std::size_t bond_index = 0; bond_index < bond_count; ++bond_index) {
        bond_centers.push_back(weighted_bond_center(molecule, geometry, bond_index));
    }

    // atom-atom part
    for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
        const auto i = static_cast<Eigen::Index>(atom_index);

        matrix(i, i) = atom_b[atom_index];
        rhs(i) = -atom_a[atom_index];

        for (std::size_t other_atom_index = atom_index + 1; other_atom_index < atom_count;
             ++other_atom_index) {
            const auto j = static_cast<Eigen::Index>(other_atom_index);
            const auto value = interaction(k, geometry.distance(atom_index, other_atom_index));

            matrix(i, j) = value;
            matrix(j, i) = value;
        }
    }

    // atom-bond part
    for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
        const auto i = static_cast<Eigen::Index>(atom_index);

        for (std::size_t bond_index = 0; bond_index < bond_count; ++bond_index) {
            const auto& bond = molecule.bond(bond_index);
            const auto bond_column = n + static_cast<Eigen::Index>(bond_index);

            if (contains_atom(bond, atom_index)) {
                matrix(i, bond_column) = atom_c[atom_index];
            } else {
                const auto value = interaction(
                    k, atom_bond_distance(geometry, atom_index, bond_centers[bond_index]));
                matrix(i, bond_column) = value;
            }
        }
    }

    // bond-atom part
    for (std::size_t bond_index = 0; bond_index < bond_count; ++bond_index) {
        const auto& bond = molecule.bond(bond_index);
        const auto row = n + static_cast<Eigen::Index>(bond_index);
        const auto& first_atom = molecule.atom(bond.first_atom_index());
        const auto& second_atom = molecule.atom(bond.second_atom_index());
        const auto first_is_canonical =
            first_atom.atomic_number() > second_atom.atomic_number() ||
            (first_atom.atomic_number() == second_atom.atomic_number() &&
             bond.first_atom_index() < bond.second_atom_index());
        const auto canonical_atom_index =
            first_is_canonical ? bond.first_atom_index() : bond.second_atom_index();

        rhs(row) = -bond_a[bond_index];

        for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
            const auto column = static_cast<Eigen::Index>(atom_index);

            if (contains_atom(bond, atom_index)) {
                if (canonical_atom_index == atom_index) {
                    matrix(row, column) = bond_d[bond_index];
                } else {
                    matrix(row, column) = bond_c[bond_index];
                }
            } else {
                const auto value = interaction(
                    k, atom_bond_distance(geometry, atom_index, bond_centers[bond_index]));
                matrix(row, column) = value;
            }
        }
    }

    // bond-bond part
    for (std::size_t bond_index = 0; bond_index < bond_count; ++bond_index) {
        const auto i = n + static_cast<Eigen::Index>(bond_index);

        matrix(i, i) = bond_b[bond_index];

        for (std::size_t other_bond_index = bond_index + 1; other_bond_index < bond_count;
             ++other_bond_index) {
            const auto j = n + static_cast<Eigen::Index>(other_bond_index);
            const auto value = interaction(
                k, core::distance(bond_centers[bond_index], bond_centers[other_bond_index]));

            matrix(i, j) = value;
            matrix(j, i) = value;
        }
    }

    Eigen::VectorXd q = solve_charge_constrained(matrix, rhs, input.target_charge());

    auto values = std::vector<double>(atom_count);

    for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
        values[atom_index] = q(static_cast<Eigen::Index>(atom_index));
    }

    for (std::size_t bond_index = 0; bond_index < bond_count; ++bond_index) {
        const auto& bond = molecule.bond(bond_index);
        const auto bond_charge = q(n + static_cast<Eigen::Index>(bond_index));

        values[bond.first_atom_index()] += 0.5 * bond_charge;
        values[bond.second_atom_index()] += 0.5 * bond_charge;
    }

    return charges::AtomicCharges{std::move(values)};
}

auto ABEEMMethod::estimate(const TargetSize& target, const MethodOptions& /*options*/) const
    -> ResourceEstimate {
    const auto size = target.atoms + target.bonds + 1.0;
    // Dense atom-plus-bond system with a charge constraint.
    return {.memory_bytes = 8.0 * size * size, .operations = 2.0 / 3.0 * size * size * size};
}

} // namespace chargefw::methods::builtin
