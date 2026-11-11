#pragma once

#include <chargefw/core/molecule.h>
#include <chargefw/parameters/models/parameter_set.h>

#include <cstddef>
#include <span>
#include <vector>

namespace chargefw::parameters {

class AtomParameterClassification {
  public:
    AtomParameterClassification() = default;

    explicit AtomParameterClassification(std::vector<std::size_t> parameter_entry_indices);

    [[nodiscard]] auto parameter_entry_indices() const noexcept -> std::span<const std::size_t>;

    [[nodiscard]] auto size() const noexcept -> std::size_t;
    [[nodiscard]] auto empty() const noexcept -> bool;

    [[nodiscard]] auto operator[](std::size_t atom_index) const noexcept -> std::size_t;
    [[nodiscard]] auto at(std::size_t atom_index) const -> std::size_t;

  private:
    std::vector<std::size_t> parameter_entry_indices_;
};

class BondParameterClassification {
  public:
    BondParameterClassification() = default;

    explicit BondParameterClassification(std::vector<std::size_t> parameter_entry_indices);

    [[nodiscard]] auto parameter_entry_indices() const noexcept -> std::span<const std::size_t>;

    [[nodiscard]] auto size() const noexcept -> std::size_t;
    [[nodiscard]] auto empty() const noexcept -> bool;

    [[nodiscard]] auto operator[](std::size_t bond_index) const noexcept -> std::size_t;
    [[nodiscard]] auto at(std::size_t bond_index) const -> std::size_t;

  private:
    std::vector<std::size_t> parameter_entry_indices_;
};

class ParameterClassification {
  public:
    ParameterClassification() = default;

    explicit ParameterClassification(AtomParameterClassification atom,
                                     BondParameterClassification bond = {});

    [[nodiscard]] auto atom() const noexcept -> const AtomParameterClassification&;
    [[nodiscard]] auto bond() const noexcept -> const BondParameterClassification&;

  private:
    AtomParameterClassification atom_;
    BondParameterClassification bond_;
};

auto validate_parameter_classification(const core::Molecule& molecule,
                                       const ParameterSet& parameters,
                                       const ParameterClassification& classification) -> void;

} // namespace chargefw::parameters
