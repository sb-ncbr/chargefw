#pragma once

#include <chargefw/io/molecule_record.h>
#include <chargefw/io/structure_options.h>

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

[[nodiscard]] auto make_record(const ::gemmi::Structure& structure,
                               std::span<const selection::SelectedModel> selected_models,
                               MoleculeRecordIdentity identity, RecordSelection selection,
                               BondStrategy bond_strategy, ConformerSelection conformer_selection,
                               std::vector<core::Bond> explicit_bonds, std::string name,
                               MolecularSourceFormat format,
                               std::vector<SourceModelMapping> source_models,
                               SourceConnectivity source_connectivity) -> ImportedMoleculeRecord;

} // namespace chargefw::io::gemmi::structure_import
