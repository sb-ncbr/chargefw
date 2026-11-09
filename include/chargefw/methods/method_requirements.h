#pragma once

#include <cstddef>
#include <cstdint>
#include <string_view>
#include <vector>

namespace chargefw::methods {

enum class ComplexityTerm : std::uint8_t {
    constant,                         // 1
    atoms,                            // N
    bonds,                            // M
    atoms_plus_bonds,                 // N + M
    atoms_squared,                    // N^2
    atoms_cubed,                      // N^3
    bonds_squared,                    // M^2
    bonds_cubed,                      // M^3
    atoms_plus_bonds_squared,         // (N + M)^2
    atoms_plus_bonds_cubed,           // (N + M)^3
    atoms_squared_plus_bonds_squared, // N^2 + M^2
};

[[nodiscard]] auto complexity_notation(ComplexityTerm value) -> std::string_view;

enum class ReducedChargePolicy : std::uint8_t {
    unsupported,
    uniform_target_global,
    zero_components,
    formal_charge_components,
    parameterized_charge_components,
};

// Size of one calculation target: a whole molecule, one conformer, or one spatial fragment.
struct TargetSize {
    double atoms = 0.0;
    double bonds = 0.0;
};

// Approximate cost of a calculation: one target in full mode for a method, or a whole plan.
// memory_bytes covers working storage beyond inputs and output charges. operations counts dense
// linear algebra in floating-point operations and other loops per atom, bond, or atom-pair visit.
struct ResourceEstimate {
    double memory_bytes = 0.0;
    double operations = 0.0;
};

struct ResourceRequirements {
    ComplexityTerm time = ComplexityTerm::constant;
    ComplexityTerm memory = ComplexityTerm::constant;

    // Reduced execution is spatial: a method may set either capability only when it also requires
    // coordinates and has a tested executor for the corresponding radius-based approximation.
    bool supports_cutoff = false;
    bool supports_cover = false;
    ReducedChargePolicy reduced_charge_policy = ReducedChargePolicy::unsupported;
};

struct MethodRequirements {
    bool coordinates = false;
    // The method uses the bonds of the molecule, not only its atoms and coordinates.
    bool bonds = false;

    std::vector<std::string_view> common_parameters;
    std::vector<std::string_view> atom_parameters;
    std::vector<std::string_view> bond_parameters;

    ResourceRequirements resources{};
    // Numerical/scientific support, independent of cutoff and cover execution support.
    bool supports_fixed_point_sources = false;

    [[nodiscard]] auto requires_common_parameters() const noexcept -> bool {
        return !common_parameters.empty();
    }

    [[nodiscard]] auto requires_atom_parameters() const noexcept -> bool {
        return !atom_parameters.empty();
    }

    [[nodiscard]] auto requires_bond_parameters() const noexcept -> bool {
        return !bond_parameters.empty();
    }

    // Atom and bond parameters are matched to each molecule; common parameters need no matching.
    [[nodiscard]] auto requires_parameter_classification() const noexcept -> bool {
        return requires_atom_parameters() || requires_bond_parameters();
    }

    [[nodiscard]] auto requires_parameters() const noexcept -> bool {
        return requires_common_parameters() || requires_atom_parameters() ||
               requires_bond_parameters();
    }
};

} // namespace chargefw::methods
