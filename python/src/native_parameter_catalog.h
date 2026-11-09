#pragma once

#include <chargefw/parameters/models/parameter_set.h>

#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace chargefw::python {

class NativeParameterCatalog {
  public:
    explicit NativeParameterCatalog(std::vector<parameters::ParameterSet> parameter_sets)
        : parameter_sets_{std::move(parameter_sets)} {}

    [[nodiscard]] auto parameter_sets() const noexcept
        -> const std::vector<parameters::ParameterSet>& {
        return parameter_sets_;
    }

    // The parameter set with the given ID; throws std::invalid_argument when there is none.
    [[nodiscard]] auto find(const std::string& id) const -> const parameters::ParameterSet& {
        for (const auto& parameter_set : parameter_sets_) {
            if (parameter_set.id() == id) {
                return parameter_set;
            }
        }
        throw std::invalid_argument{"no parameter set '" + id + "' in the catalog"};
    }

  private:
    std::vector<parameters::ParameterSet> parameter_sets_;
};

} // namespace chargefw::python
