#pragma once

#include <chargefw/core/position.h>
#include <chargefw/features/conformer_features.h>
#include <chargefw/features/prepared_molecule.h>
#include <chargefw/methods/method_options.h>
#include <chargefw/parameters/models/parameter_view.h>

#include <functional>
#include <span>
#include <stdexcept>

namespace chargefw::methods {

struct FixedPointSource {
    core::Position position;
    double charge = 0.0;
};

class CalculationInput {
  public:
    // fixed_sources and reference_charges are borrowed and must outlive every call using this
    // input.
    //
    // reference_charges is empty or has one value per atom. A method that starts from per-atom
    // reference charges uses these values as they are, instead of deriving and normalizing its
    // own, and its result keeps their sum instead of target_charge. Methods without reference
    // charges ignore them. Reduced execution passes each fragment the values of its atoms
    // computed for the whole molecule; they need not sum to the fragment's formal charge.
    CalculationInput(const features::PreparedMolecule& prepared_molecule,
                     const MethodOptions& method_options, double target_charge,
                     const features::ConformerFeatures* geometry = nullptr,
                     const parameters::ParameterView* parameters = nullptr,
                     std::span<const FixedPointSource> fixed_sources = {},
                     std::span<const double> reference_charges = {});

    CalculationInput(features::PreparedMolecule&&, const MethodOptions&, double,
                     const features::ConformerFeatures* = nullptr,
                     const parameters::ParameterView* = nullptr,
                     std::span<const FixedPointSource> = {}, std::span<const double> = {}) = delete;

    CalculationInput(const features::PreparedMolecule&, MethodOptions&&, double,
                     const features::ConformerFeatures* = nullptr,
                     const parameters::ParameterView* = nullptr,
                     std::span<const FixedPointSource> = {}, std::span<const double> = {}) = delete;

    [[nodiscard]] auto prepared_molecule() const noexcept -> const features::PreparedMolecule&;

    [[nodiscard]] auto molecule() const noexcept -> const core::Molecule&;

    [[nodiscard]] auto topology() const noexcept -> const features::TopologyFeatures&;

    [[nodiscard]] auto method_options() const noexcept -> const MethodOptions&;

    [[nodiscard]] auto target_charge() const noexcept -> double;

    [[nodiscard]] auto fixed_sources() const noexcept -> std::span<const FixedPointSource>;

    [[nodiscard]] auto reference_charges() const noexcept -> std::span<const double>;

    [[nodiscard]] auto geometry() const -> const features::ConformerFeatures&;

    [[nodiscard]] auto parameters() const -> const parameters::ParameterView&;

  private:
    std::reference_wrapper<const features::PreparedMolecule> prepared_molecule_;
    std::reference_wrapper<const MethodOptions> method_options_;
    double target_charge_;
    const features::ConformerFeatures* geometry_ = nullptr;
    const parameters::ParameterView* parameters_ = nullptr;
    std::span<const FixedPointSource> fixed_sources_;
    std::span<const double> reference_charges_;
};

} // namespace chargefw::methods
