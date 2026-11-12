#pragma once

#include <chargefw/parameters/classification/parameter_classification.h>
#include <chargefw/parameters/models/parameter_set.h>

#include <cstddef>
#include <string>
#include <string_view>

namespace chargefw::parameters {

class ParameterView;

class AtomParameterAccessor {
  public:
    [[nodiscard]] auto operator[](std::size_t atom_index) const -> double;
    [[nodiscard]] auto at(std::size_t atom_index) const -> double;

  private:
    friend class ParameterView;

    AtomParameterAccessor(const AtomParameters& parameters,
                          const AtomParameterClassification& classification,
                          std::string parameter_name);

    const AtomParameters* parameters_;
    const AtomParameterClassification* classification_;
    std::string parameter_name_;
};

class BondParameterAccessor {
  public:
    [[nodiscard]] auto operator[](std::size_t bond_index) const -> double;
    [[nodiscard]] auto at(std::size_t bond_index) const -> double;

  private:
    friend class ParameterView;

    BondParameterAccessor(const BondParameters& parameters,
                          const BondParameterClassification& classification,
                          std::string parameter_name);

    const BondParameters* parameters_;
    const BondParameterClassification* classification_;
    std::string parameter_name_;
};

class ParameterView {
  public:
    explicit ParameterView(const ParameterSet& parameter_set,
                           const ParameterClassification& classification);

    explicit ParameterView(const ParameterSet&& parameter_set,
                           const ParameterClassification& classification) = delete;
    explicit ParameterView(const ParameterSet& parameter_set,
                           const ParameterClassification&& classification) = delete;
    explicit ParameterView(const ParameterSet&& parameter_set,
                           const ParameterClassification&& classification) = delete;

    [[nodiscard]] auto parameter_set() const noexcept -> const ParameterSet&;
    [[nodiscard]] auto classification() const noexcept -> const ParameterClassification&;

    [[nodiscard]] auto common(std::string_view name) const -> double;

    [[nodiscard]] auto atom(std::string_view name) const -> AtomParameterAccessor;
    [[nodiscard]] auto bond(std::string_view name) const -> BondParameterAccessor;

  private:
    const ParameterSet* parameter_set_;
    const ParameterClassification* classification_;
};

} // namespace chargefw::parameters
