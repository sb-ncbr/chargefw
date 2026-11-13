#pragma once

#include <chargefw/io/molecule_record.h>

#include <optional>
#include <utility>
#include <vector>

namespace chargefw::python {

class NativeInputMetadata {
  public:
    NativeInputMetadata(io::MoleculeRecordIdentity identity, std::vector<io::ImportIssue> issues,
                        std::optional<io::MoleculeImportMetadata> import_metadata)
        : identity_{std::move(identity)}, issues_{std::move(issues)},
          import_metadata_{std::move(import_metadata)} {}

    [[nodiscard]] auto make_record(core::Molecule molecule) const -> io::ImportedMoleculeRecord {
        return io::ImportedMoleculeRecord{.molecule = std::move(molecule),
                                          .identity = identity_,
                                          .caller_atom_ids = std::nullopt,
                                          .import_warnings = issues_,
                                          .import_metadata = import_metadata_};
    }

  private:
    io::MoleculeRecordIdentity identity_;
    std::vector<io::ImportIssue> issues_;
    std::optional<io::MoleculeImportMetadata> import_metadata_;
};

} // namespace chargefw::python
