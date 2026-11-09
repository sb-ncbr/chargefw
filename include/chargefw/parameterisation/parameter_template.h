#pragma once

#include <chargefw/features/prepared_molecule_collection.h>
#include <chargefw/parameters/models/parameter_key.h>
#include <chargefw/parameters/models/parameter_set.h>
#include <chargefw/parameters/models/parameter_set_metadata.h>

#include <string>
#include <vector>

namespace chargefw::parameterisation {

// The parameter names a template gives its common parameters and each of its atom and bond entries.
struct ParameterNames {
    std::vector<std::string> common;
    std::vector<std::string> atom;
    std::vector<std::string> bond;
};

// A parameter set with an entry for every atom type and bond type that occurs in `molecules` under
// the given classifiers, the given parameter names, and all values zero. Atom entries are made only
// when `names.atom` is not empty, and bond entries only when `names.bond` is not empty; the atom
// keys of a bond entry use `atom_types`. Entries are ordered by atomic number, then type, so the
// same molecules and classifiers always give the same layout. Every atom and bond of `molecules`
// matches an entry under strict matching. `metadata` must name the method.
[[nodiscard]] auto make_parameter_template(const features::PreparedMoleculeCollection& molecules,
                                           const ParameterNames& names,
                                           parameters::AtomParameterClassificationKind atom_types,
                                           parameters::BondParameterClassificationKind bond_types,
                                           parameters::ParameterSetMetadata metadata)
    -> parameters::ParameterSet;

} // namespace chargefw::parameterisation
