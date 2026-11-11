#include <chargefw/parameters/classification/parameter_classification.h>

#include <stdexcept>
#include <string>
#include <utility>

namespace chargefw::parameters {
namespace {

auto validate_mapping_size(std::span<const std::size_t> indices, const std::size_t expected_size,
                           const bool mapping_required, const std::string_view context) -> void {
    if (!mapping_required) {
        if (!indices.empty()) {
            throw std::invalid_argument{
                std::string{context} + " classification is present but no parameter entries exist"};
        }
        return;
    }

    if (indices.size() != expected_size) {
        throw std::invalid_argument{std::string{context} +
                                    " classification size does not match molecule"};
    }
}

auto validate_parameter_entry_indices(std::span<const std::size_t> indices,
                                      const std::size_t parameter_entry_count,
                                      const std::string_view context) -> void {
    for (const auto parameter_entry_index : indices) {
        if (parameter_entry_index >= parameter_entry_count) {
            throw std::invalid_argument{std::string{context} +
                                        " classification refers to an unknown parameter entry"};
        }
    }
}

} // namespace

AtomParameterClassification::AtomParameterClassification(
    std::vector<std::size_t> parameter_entry_indices)
    : parameter_entry_indices_{std::move(parameter_entry_indices)} {}

auto AtomParameterClassification::parameter_entry_indices() const noexcept
    -> std::span<const std::size_t> {
    return parameter_entry_indices_;
}

auto AtomParameterClassification::size() const noexcept -> std::size_t {
    return parameter_entry_indices_.size();
}

auto AtomParameterClassification::empty() const noexcept -> bool {
    return parameter_entry_indices_.empty();
}

auto AtomParameterClassification::operator[](const std::size_t atom_index) const noexcept
    -> std::size_t {
    return parameter_entry_indices_[atom_index];
}

auto AtomParameterClassification::at(const std::size_t atom_index) const -> std::size_t {
    return parameter_entry_indices_.at(atom_index);
}

BondParameterClassification::BondParameterClassification(
    std::vector<std::size_t> parameter_entry_indices)
    : parameter_entry_indices_{std::move(parameter_entry_indices)} {}

auto BondParameterClassification::parameter_entry_indices() const noexcept
    -> std::span<const std::size_t> {
    return parameter_entry_indices_;
}

auto BondParameterClassification::size() const noexcept -> std::size_t {
    return parameter_entry_indices_.size();
}

auto BondParameterClassification::empty() const noexcept -> bool {
    return parameter_entry_indices_.empty();
}

auto BondParameterClassification::operator[](const std::size_t bond_index) const noexcept
    -> std::size_t {
    return parameter_entry_indices_[bond_index];
}

auto BondParameterClassification::at(const std::size_t bond_index) const -> std::size_t {
    return parameter_entry_indices_.at(bond_index);
}

ParameterClassification::ParameterClassification(AtomParameterClassification atom,
                                                 BondParameterClassification bond)
    : atom_{std::move(atom)}, bond_{std::move(bond)} {}

auto ParameterClassification::atom() const noexcept -> const AtomParameterClassification& {
    return atom_;
}

auto ParameterClassification::bond() const noexcept -> const BondParameterClassification& {
    return bond_;
}

auto validate_parameter_classification(const core::Molecule& molecule,
                                       const ParameterSet& parameters,
                                       const ParameterClassification& classification) -> void {
    const auto atom_indices = classification.atom().parameter_entry_indices();

    validate_mapping_size(atom_indices, molecule.atom_count(), !parameters.atom().empty(), "atom");
    validate_parameter_entry_indices(atom_indices, parameters.atom().size(), "atom");

    const auto bond_indices = classification.bond().parameter_entry_indices();

    validate_mapping_size(bond_indices, molecule.bond_count(), !parameters.bond().empty(), "bond");
    validate_parameter_entry_indices(bond_indices, parameters.bond().size(), "bond");
}

} // namespace chargefw::parameters
