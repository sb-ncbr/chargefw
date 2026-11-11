#include "methods/builtin/eem.h"
#include "methods/builtin/dense_solve.h"
#include "methods/builtin/eem_system.h"
#include "methods/fixed_source_validation.h"

#include <chargefw/core/molecule.h>
#include <chargefw/core/position.h>
#include <chargefw/parameters/models/parameter_view.h>

#include <Eigen/Core>

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace chargefw::methods::builtin {
namespace {

[[nodiscard]] auto fixed_source_potential(const CalculationInput& input, const double kappa)
    -> std::vector<double> {
    const auto sources = input.fixed_sources();
    if (sources.empty()) {
        return {};
    }

    const auto& geometry = input.geometry();
    auto potential = std::vector<double>(input.molecule().atom_count(), 0.0);
    for (std::size_t atom_index = 0; atom_index < potential.size(); ++atom_index) {
        const auto& active_position = geometry.position(atom_index);
        for (const auto& source : sources) {
            const auto distance = core::distance(active_position, source.position);
            potential[atom_index] += kappa * source.charge / distance;
        }
    }
    return potential;
}

} // namespace

auto eem_core::build_system(const CalculationInput& input) -> System {
    const auto& molecule = input.molecule();
    const auto& parameters = input.parameters();
    const auto& geometry = input.geometry();
    const auto atom_count = molecule.atom_count();

    const auto kappa = parameters.common("kappa");
    if (kappa <= 0.0) {
        throw std::logic_error{"EEM common parameter 'kappa' must be positive"};
    }
    const auto source_potential = fixed_source_potential(input, kappa);

    const auto A_parameter = parameters.atom("A");
    const auto B_parameter = parameters.atom("B");

    const auto n = static_cast<Eigen::Index>(atom_count);

    auto system =
        System{.matrix = Eigen::MatrixXd::Zero(n + 1, n + 1), .rhs = Eigen::VectorXd::Zero(n + 1)};
    auto& A = system.matrix;
    auto& b = system.rhs;

    for (std::size_t atom_index = 0; atom_index < atom_count; ++atom_index) {
        const auto i = static_cast<Eigen::Index>(atom_index);

        A(i, i) = B_parameter[atom_index];
        b(i) = -A_parameter[atom_index];
        if (!source_potential.empty()) {
            b(i) -= source_potential[atom_index];
        }

        for (std::size_t other_atom_index = atom_index + 1; other_atom_index < atom_count;
             ++other_atom_index) {
            const auto j = static_cast<Eigen::Index>(other_atom_index);
            const auto distance = geometry.distance(atom_index, other_atom_index);
            const auto value = kappa / distance;

            A(i, j) = value;
            A(j, i) = value;
        }
    }

    set_charge_constraint(A, b, input.target_charge());

    return system;
}

auto EEMMethod::calculate(const CalculationInput& input) const -> charges::AtomicCharges {
    detail::validate_fixed_source_input(input, "EEM");
    const auto atom_count = input.molecule().atom_count();

    if (atom_count == 0) {
        return charges::AtomicCharges{std::vector<double>{}};
    }

    auto system = eem_core::build_system(input);
    const Eigen::VectorXd q =
        solve_in_place(system.matrix, system.rhs).head(static_cast<Eigen::Index>(atom_count));

    return charges::AtomicCharges{std::vector<double>{q.data(), q.data() + q.size()}};
}

auto EEMMethod::estimate(const TargetSize& target, const MethodOptions& /*options*/) const
    -> ResourceEstimate {
    return charge_constrained_estimate(target.atoms);
}

} // namespace chargefw::methods::builtin
