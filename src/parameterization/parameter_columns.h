#pragma once

#include "parameterization/method_derivatives.h"

#include <Eigen/Core>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace chargefw::parameterization::detail {

inline auto require_finite(const std::span<const double> values, const std::string_view message)
    -> void {
    if (!std::ranges::all_of(values, [](const double value) { return std::isfinite(value); })) {
        throw std::runtime_error{std::string{message}};
    }
}

[[nodiscard]] inline auto valid_column_assignment(const std::span<const std::size_t> assignment,
                                                  const std::size_t size,
                                                  const std::size_t column_count) -> bool {
    return assignment.size() == size &&
           std::ranges::all_of(assignment, [column_count](const std::size_t column) {
               return column < column_count || column == ParameterColumns::none;
           });
}

// Throws std::invalid_argument unless `columns` assigns every value of a method with these
// parameters, on a molecule with these atoms and bonds, to a column below columns.count or to none.
inline auto check_columns(const ParameterColumns& columns, const std::size_t atom_parameters,
                          const std::size_t bond_parameters, const std::size_t common_parameters,
                          const std::size_t atom_count, const std::size_t bond_count,
                          const std::string& method) -> void {
    const auto valid = [&columns](const std::span<const std::size_t> assignment,
                                  const std::size_t size) {
        return valid_column_assignment(assignment, size, columns.count);
    };
    const auto all_valid = [&valid](const std::vector<std::vector<std::size_t>>& channels,
                                    const std::size_t size) {
        return std::ranges::all_of(channels,
                                   [&valid, size](const std::vector<std::size_t>& channel) {
                                       return valid(channel, size);
                                   });
    };
    if (columns.atom.size() != atom_parameters || columns.bond.size() != bond_parameters ||
        !all_valid(columns.atom, atom_count) || !all_valid(columns.bond, bond_count) ||
        !valid(columns.common, common_parameters)) {
        throw std::invalid_argument{method +
                                    " jacobian: the columns do not match the method's parameters"};
    }
}

// Some columns of a jacobian() as a dense block: block column j belongs to columns[j].
struct ColumnBlock {
    std::vector<std::size_t> columns;
    Eigen::MatrixXd values;
};

// The columns that the values of one parameter belong to, ascending, with value i's row of the
// block 1 in its column and 0 elsewhere.
[[nodiscard]] inline auto indicator_block(const std::span<const std::size_t> assignment)
    -> ColumnBlock {
    auto block = ColumnBlock{};
    for (const auto column : assignment) {
        if (column != ParameterColumns::none) {
            block.columns.push_back(column);
        }
    }
    std::ranges::sort(block.columns);
    const auto [first, last] = std::ranges::unique(block.columns);
    block.columns.erase(first, last);

    block.values = Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(assignment.size()),
                                         static_cast<Eigen::Index>(block.columns.size()));
    for (std::size_t index = 0; index < assignment.size(); ++index) {
        if (assignment[index] != ParameterColumns::none) {
            const auto position = std::ranges::lower_bound(block.columns, assignment[index]);
            block.values(static_cast<Eigen::Index>(index),
                         static_cast<Eigen::Index>(position - block.columns.begin())) = 1.0;
        }
    }
    return block;
}

} // namespace chargefw::parameterization::detail
