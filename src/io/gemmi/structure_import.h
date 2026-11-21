#pragma once

#include <chargefw/io/molecule_record.h>
#include <chargefw/io/structure_options.h>

#include "bonds.h"
#include "selection.h"

#include <gemmi/model.hpp>

#include <cstddef>
#include <span>
#include <string>
#include <vector>

namespace chargefw::io::gemmi::structure_import {

struct SourceModelMapping {
    SourceConformerReference conformer;
};

// explicit_bonds are the bonds read from the file, when the bond strategy reads them. The record
// receives the warnings about them too.
[[nodiscard]] auto
make_record(const ::gemmi::Structure& structure,
            std::span<const selection::SelectedModel> selected_models,
            MoleculeRecordIdentity identity, RecordSelection selection, BondStrategy bond_strategy,
            const ComponentDictionary* dictionary, ConformerSelection conformer_selection,
            bonds::ExplicitBonds explicit_bonds, std::string name, MolecularSourceFormat format,
            std::vector<SourceModelMapping> source_models, SourceConnectivity source_connectivity)
    -> ImportedMoleculeRecord;

} // namespace chargefw::io::gemmi::structure_import
