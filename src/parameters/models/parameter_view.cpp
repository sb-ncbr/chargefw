#include <chargefw/parameters/models/parameter_view.h>

#include <utility>

namespace chargefw::parameters {

AtomParameterAccessor::AtomParameterAccessor(const AtomParameters& parameters,
                                             const AtomParameterClassification& classification,
                                             std::string parameter_name)
    : parameters_{&parameters}, classification_{&classification},
      parameter_name_{std::move(parameter_name)} {}

auto AtomParameterAccessor::operator[](const std::size_t atom_index) const -> double {
    const auto parameter_entry_index = (*classification_)[atom_index];
    return parameters_->parameter(parameter_entry_index, parameter_name_);
}

auto AtomParameterAccessor::at(const std::size_t atom_index) const -> double {
    const auto parameter_entry_index = classification_->at(atom_index);
    return parameters_->parameter(parameter_entry_index, parameter_name_);
}

BondParameterAccessor::BondParameterAccessor(const BondParameters& parameters,
                                             const BondParameterClassification& classification,
                                             std::string parameter_name)
    : parameters_{&parameters}, classification_{&classification},
      parameter_name_{std::move(parameter_name)} {}

auto BondParameterAccessor::operator[](const std::size_t bond_index) const -> double {
    const auto parameter_entry_index = (*classification_)[bond_index];
    return parameters_->parameter(parameter_entry_index, parameter_name_);
}

auto BondParameterAccessor::at(const std::size_t bond_index) const -> double {
    const auto parameter_entry_index = classification_->at(bond_index);
    return parameters_->parameter(parameter_entry_index, parameter_name_);
}

ParameterView::ParameterView(const ParameterSet& parameter_set,
                             const ParameterClassification& classification)
    : parameter_set_{&parameter_set}, classification_{&classification} {}

auto ParameterView::parameter_set() const noexcept -> const ParameterSet& {
    return *parameter_set_;
}

auto ParameterView::classification() const noexcept -> const ParameterClassification& {
    return *classification_;
}

auto ParameterView::common(const std::string_view name) const -> double {
    return parameter_set_->common().parameter(name);
}

auto ParameterView::atom(const std::string_view name) const -> AtomParameterAccessor {
    return AtomParameterAccessor{parameter_set_->atom(), classification_->atom(),
                                 std::string{name}};
}

auto ParameterView::bond(const std::string_view name) const -> BondParameterAccessor {
    return BondParameterAccessor{parameter_set_->bond(), classification_->bond(),
                                 std::string{name}};
}

} // namespace chargefw::parameters
