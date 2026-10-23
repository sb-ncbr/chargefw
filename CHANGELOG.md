# Changelog

## 3.0.0a1

First pre-release of the ChargeFW 3 series, the successor of
[ChargeFW2](https://github.com/sb-ncbr/ChargeFW2). The version number continues ChargeFW2's lineage;
the `chargefw` Python package moves from 0.1.7 directly to 3.0.0a1. Install it with
`pip install --pre chargefw`.

ChargeFW 3 is a new implementation. Its highlights are:

- one native engine behind a C++23 library, a command-line application, and a Python package for
  CPython 3.11 to 3.14 on Linux x86-64;
- 22 methods and the bundled parameter sets, with deterministic method and parameter-set selection and
  structured reasons for rejected candidates;
- explicit molecular input: topology, protonation, formal charges, and coordinates are never changed
  silently, and missing prerequisites are reported as failures;
- full calculations and explicit cutoff and cover approximations, recorded in result provenance;
- SDF, MOL2, PDB, mmCIF, and molecule JSON input; result JSON, MOL2, and mmCIF output; source atom,
  molecule, and conformer order preserved in every result;
- optional Gemmi and RDKit integration in Python.

Documentation is published at [chargefw.readthedocs.io](https://chargefw.readthedocs.io/).

### Changes since 0.1.7

- Parameter-set keys are exposed in Python (`atom_keys` and `bond_keys`) and in the CLI
  (`chargefw parameters ID --keys`).
- The bundled periodic table is exposed in Python as `chargefw.periodic_table`.
- New recipes calculate dipole moments and attach charges to RDKit SDF records; the file calculation
  recipe can export CSV, and the fixed-ion and molecule-inspection recipes report more detail.
- The repository moved to [sb-ncbr/chargefw](https://github.com/sb-ncbr/chargefw).
