#pragma once

#include <chargefw/parameters/models/atom_parameters.h>
#include <chargefw/parameters/models/bond_parameters.h>
#include <chargefw/parameters/models/common_parameters.h>
#include <chargefw/parameters/models/parameter_set_metadata.h>

#include <string_view>

namespace chargefw::parameters {

class ParameterSet {
  public:
    explicit ParameterSet(ParameterSetMetadata metadata, CommonParameters common = {},
                          AtomParameters atom = {}, BondParameters bond = {});

    [[nodiscard]] auto metadata() const noexcept -> const ParameterSetMetadata&;

    [[nodiscard]] auto id() const noexcept -> std::string_view;
    [[nodiscard]] auto method_id() const noexcept -> std::string_view;
    [[nodiscard]] auto name() const noexcept -> std::string_view;
    [[nodiscard]] auto publication() const noexcept -> std::string_view;
    [[nodiscard]] auto notes() const noexcept -> std::string_view;
    [[nodiscard]] auto priority() const noexcept -> std::uint16_t;

    [[nodiscard]] auto common() const noexcept -> const CommonParameters&;
    [[nodiscard]] auto atom() const noexcept -> const AtomParameters&;
    [[nodiscard]] auto bond() const noexcept -> const BondParameters&;

  private:
    ParameterSetMetadata metadata_;
    CommonParameters common_;
    AtomParameters atom_;
    BondParameters bond_;
};

// True when `candidate` has the method ID, entry keys, and parameter names of `reference`, in the
// same order; the values and the descriptive metadata may differ.
[[nodiscard]] auto compatible(const ParameterSet& reference, const ParameterSet& candidate) -> bool;

} // namespace chargefw::parameters