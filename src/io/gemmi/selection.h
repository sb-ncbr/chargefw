#pragma once

#include <chargefw/io/structure_options.h>

#include "io/component.h"

#include <gemmi/model.hpp>

#include <cstddef>
#include <functional>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace chargefw::io::gemmi::selection {

struct SelectedResidue {
    const ::gemmi::Residue* residue;
    std::string_view chain_name;
    std::string component_name;
    std::vector<detail::ComponentAtom> atoms;

    [[nodiscard]] auto component() const -> detail::ComponentView {
        return {.name = component_name, .atoms = atoms};
    }
};

// The name that tells the atoms of a residue apart, where it is not the one Gemmi gives the atom.
using AtomName = std::function<std::string(const ::gemmi::Atom&)>;

// Borrows a Gemmi model and its atoms. The model must outlive this view and remain unmodified.
class SelectedModel {
  public:
    SelectedModel(const ::gemmi::Model& model, RecordSelection selection,
                  const AtomName& atom_name = {});

    [[nodiscard]] auto atoms() const -> const std::vector<const ::gemmi::Atom*>&;
    [[nodiscard]] auto sites() const -> const std::vector<::gemmi::const_CRA>&;
    [[nodiscard]] auto residues() const -> const std::vector<SelectedResidue>&;
    [[nodiscard]] auto atom_index_by_serial(int serial) const -> std::optional<std::size_t>;
    void normalize_component_names(std::span<const SourceAtomReference> source_atoms);

  private:
    std::vector<const ::gemmi::Atom*> atoms_;
    std::vector<::gemmi::const_CRA> sites_;
    std::vector<SelectedResidue> residues_;
    std::unordered_map<int, std::size_t> serial_indices_;
};

[[nodiscard]] auto select_models(const ::gemmi::Structure& structure, RecordSelection selection,
                                 ConformerSelection conformers, const AtomName& atom_name = {})
    -> std::vector<SelectedModel>;

} // namespace chargefw::io::gemmi::selection
