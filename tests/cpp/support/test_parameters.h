#pragma once

#include <chargefw/parameters/models/atom_parameters.h>
#include <chargefw/parameters/models/bond_parameters.h>
#include <chargefw/parameters/models/common_parameters.h>
#include <chargefw/parameters/models/parameter_key.h>
#include <chargefw/parameters/models/parameter_set.h>
#include <chargefw/parameters/models/parameter_set_metadata.h>

#include <optional>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace chargefw::test {

[[nodiscard]] inline auto atom_key(const int atomic_number,
                                   const parameters::AtomParameterClassificationKind classification,
                                   std::string type) -> parameters::AtomParameterKey {
    return {
        .atomic_number = atomic_number, .classification = classification, .type = std::move(type)};
}

[[nodiscard]] inline auto plain_atom_key(const int atomic_number) -> parameters::AtomParameterKey {
    return atom_key(atomic_number, parameters::AtomParameterClassificationKind::PLAIN, "*");
}

[[nodiscard]] inline auto hbo_atom_key(const int atomic_number, std::string type)
    -> parameters::AtomParameterKey {
    return atom_key(atomic_number, parameters::AtomParameterClassificationKind::HIGHEST_BOND_ORDER,
                    std::move(type));
}

[[nodiscard]] inline auto single_bond_key(const int first_atomic_number,
                                          const int second_atomic_number)
    -> parameters::BondParameterKey {
    return {.first_atom = plain_atom_key(first_atomic_number),
            .second_atom = plain_atom_key(second_atomic_number),
            .bond = {.classification = parameters::BondParameterClassificationKind::BOND_ORDER,
                     .type = "1"}};
}

/// Undirected plain bond key: matches a bond between the two elements regardless of direction.
[[nodiscard]] inline auto plain_bond_key(const int first_atomic_number,
                                         const int second_atomic_number)
    -> parameters::BondParameterKey {
    return {.first_atom = plain_atom_key(first_atomic_number),
            .second_atom = plain_atom_key(second_atomic_number),
            .bond = {.classification = parameters::BondParameterClassificationKind::PLAIN,
                     .type = "*"}};
}

/// Wildcard plain bond key: matches any bond regardless of element or direction.
[[nodiscard]] inline auto plain_bond_key() -> parameters::BondParameterKey {
    return {.first_atom = plain_atom_key(0),
            .second_atom = plain_atom_key(0),
            .bond = {.classification = parameters::BondParameterClassificationKind::PLAIN,
                     .type = "*"}};
}

[[nodiscard]] inline auto make_qeq_ho_parameters() -> parameters::ParameterSet {
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{
            .id = "test-qeq", .method_id = "qeq", .name = "Test QEq parameters"},
        {},
        parameters::AtomParameters{{{.key = plain_atom_key(1),
                                     .parameters = {{.name = "electronegativity", .value = 4.5280},
                                                    {.name = "hardness", .value = 13.8904}}},
                                    {.key = plain_atom_key(8),
                                     .parameters = {{.name = "electronegativity", .value = 8.741},
                                                    {.name = "hardness", .value = 13.364}}}}}};
}

[[nodiscard]] inline auto make_sfkeem_ho_parameters() -> parameters::ParameterSet {
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{
            .id = "test-sfkeem", .method_id = "sfkeem", .name = "Test SFKEEM parameters"},
        parameters::CommonParameters{{{.name = "sigma", .value = 1.0}}},
        parameters::AtomParameters{
            {{.key = plain_atom_key(1),
              .parameters = {{.name = "A", .value = 1.0}, {.name = "B", .value = 10.0}}},
             {.key = plain_atom_key(8),
              .parameters = {{.name = "A", .value = 2.0}, {.name = "B", .value = 10.0}}}}}};
}

[[nodiscard]] inline auto make_abeem_ho_parameters() -> parameters::ParameterSet {
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{
            .id = "test-abeem", .method_id = "abeem", .name = "Test ABEEM parameters"},
        parameters::CommonParameters{{{.name = "k", .value = 1.0}}},
        parameters::AtomParameters{{{.key = plain_atom_key(1),
                                     .parameters = {{.name = "a", .value = 1.0},
                                                    {.name = "b", .value = 10.0},
                                                    {.name = "c", .value = 0.5}}},
                                    {.key = plain_atom_key(8),
                                     .parameters = {{.name = "a", .value = 2.0},
                                                    {.name = "b", .value = 10.0},
                                                    {.name = "c", .value = 0.5}}}}},
        parameters::BondParameters{{{.key = plain_bond_key(8, 1),
                                     .parameters = {{.name = "A", .value = 1.0},
                                                    {.name = "B", .value = 10.0},
                                                    {.name = "C", .value = 0.5},
                                                    {.name = "D", .value = 0.5}}}}}};
}

/// EEM parameters for H and O with equal hardness.
[[nodiscard]] inline auto make_eem_parameters() -> parameters::ParameterSet {
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{.id = "test-eem", .method_id = "eem", .name = "Test EEM"},
        parameters::CommonParameters{{{.name = "kappa", .value = 1.0}}},
        parameters::AtomParameters{
            {{.key = plain_atom_key(1),
              .parameters = {{.name = "A", .value = 1.0}, {.name = "B", .value = 10.0}}},
             {.key = plain_atom_key(8),
              .parameters = {{.name = "A", .value = 2.0}, {.name = "B", .value = 10.0}}}}}};
}

/// EEM parameters for H (A = 1, B = 5) and O (A = 2, B = 9) with a selectable kappa, chosen so
/// that two-atom solutions have simple closed forms.
[[nodiscard]] inline auto make_eem_ho_parameters(const double kappa) -> parameters::ParameterSet {
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{
            .id = "test-eem-ho", .method_id = "eem", .name = "Test EEM H/O"},
        parameters::CommonParameters{{{.name = "kappa", .value = kappa}}},
        parameters::AtomParameters{
            {{.key = plain_atom_key(1),
              .parameters = {{.name = "A", .value = 1.0}, {.name = "B", .value = 5.0}}},
             {.key = plain_atom_key(8),
              .parameters = {{.name = "A", .value = 2.0}, {.name = "B", .value = 9.0}}}}}};
}

/// Per-element values that vary between SQE-family H/O fixtures.
struct SqeHoValues {
    double hydrogen_width = 1.0;
    double oxygen_width = 1.0;
    std::optional<double> hydrogen_q0;
    std::optional<double> oxygen_q0;
    double kappa = 1.0;
};

/// SQE-family (sqe, sqeq0, sqeqp) parameters for H and O with an H-O single-bond kappa of 1 by
/// default.
[[nodiscard]] inline auto make_sqe_ho_parameters(const std::string_view method_id,
                                                 const SqeHoValues& values = {})
    -> parameters::ParameterSet {
    auto hydrogen =
        std::vector<parameters::NamedParameter>{{.name = "electronegativity", .value = 4.5280},
                                                {.name = "hardness", .value = 13.8904},
                                                {.name = "width", .value = values.hydrogen_width}};
    auto oxygen =
        std::vector<parameters::NamedParameter>{{.name = "electronegativity", .value = 8.741},
                                                {.name = "hardness", .value = 13.364},
                                                {.name = "width", .value = values.oxygen_width}};
    if (values.hydrogen_q0.has_value()) {
        hydrogen.push_back({.name = "q0", .value = *values.hydrogen_q0});
    }
    if (values.oxygen_q0.has_value()) {
        oxygen.push_back({.name = "q0", .value = *values.oxygen_q0});
    }
    return parameters::ParameterSet{
        parameters::ParameterSetMetadata{.id = "test-" + std::string{method_id},
                                         .method_id = std::string{method_id},
                                         .name = "Test SQE-family parameters"},
        {},
        parameters::AtomParameters{{{.key = plain_atom_key(1), .parameters = std::move(hydrogen)},
                                    {.key = plain_atom_key(8), .parameters = std::move(oxygen)}}},
        parameters::BondParameters{{{.key = single_bond_key(1, 8),
                                     .parameters = {{.name = "kappa", .value = values.kappa}}}}}};
}

} // namespace chargefw::test
