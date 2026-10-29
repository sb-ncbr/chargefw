#pragma once

#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/conformer.h>

#include <cstddef>
#include <span>
#include <string>
#include <string_view>
#include <vector>

namespace chargefw::core {
class Molecule {
  public:
    explicit Molecule(std::vector<Atom> atoms, std::vector<Bond> bonds = {},
                      std::vector<Conformer> conformers = {}, std::string name = {});

    [[nodiscard]] auto name() const noexcept -> std::string_view;

    [[nodiscard]] auto atoms() const noexcept -> std::span<const Atom>;
    [[nodiscard]] auto bonds() const noexcept -> std::span<const Bond>;
    [[nodiscard]] auto conformers() const noexcept -> std::span<const Conformer>;

    [[nodiscard]] auto atom_count() const noexcept -> std::size_t;
    [[nodiscard]] auto bond_count() const noexcept -> std::size_t;
    [[nodiscard]] auto conformer_count() const noexcept -> std::size_t;

    [[nodiscard]] auto atom(std::size_t index) const -> const Atom&;
    [[nodiscard]] auto bond(std::size_t index) const -> const Bond&;
    [[nodiscard]] auto conformer(std::size_t index) const -> const Conformer&;

    [[nodiscard]] auto has_coordinates() const noexcept -> bool;

    [[nodiscard]] auto operator==(const Molecule&) const -> bool = default;

  private:
    std::vector<Atom> atoms_;
    std::vector<Bond> bonds_;
    std::vector<Conformer> conformers_;

    std::string name_;
};

[[nodiscard]] auto total_formal_charge(const Molecule& molecule) noexcept -> double;

} // namespace chargefw::core