#pragma once

#include <chargefw/parameters/models/common_parameters.h>
#include <chargefw/parameters/models/parameter_key.h>

#include <cstddef>
#include <span>
#include <string_view>
#include <vector>

namespace chargefw::parameters {

struct BondParameterEntry {
    BondParameterKey key;
    std::vector<NamedParameter> parameters;
};

class BondParameters {
  public:
    BondParameters() = default;

    explicit BondParameters(std::vector<BondParameterEntry> entries);

    [[nodiscard]] auto entries() const noexcept -> std::span<const BondParameterEntry>;

    [[nodiscard]] auto size() const noexcept -> std::size_t;
    [[nodiscard]] auto empty() const noexcept -> bool;

    [[nodiscard]] auto operator[](std::size_t index) const noexcept -> const BondParameterEntry&;
    [[nodiscard]] auto at(std::size_t index) const -> const BondParameterEntry&;

    [[nodiscard]] auto contains(std::size_t entry_index, std::string_view name) const -> bool;

    [[nodiscard]] auto parameter(std::size_t entry_index, std::string_view name) const -> double;

  private:
    std::vector<BondParameterEntry> entries_;
};

} // namespace chargefw::parameters
