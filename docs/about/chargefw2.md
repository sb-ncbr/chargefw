# Differences from ChargeFW2

ChargeFW is the successor of [ChargeFW2](https://github.com/sb-ncbr/ChargeFW2), the program behind
[Atomic Charge Calculator](https://acc.biodata.ceitec.cz/). It keeps the methods and parameter sets of
ChargeFW2 and its version numbers continue from it, but it is a new implementation with a different
command line and Python interface.

## What is new

- A Python package that works with molecules in memory: read files, build molecules from NumPy arrays,
  or convert [RDKit](../interfaces/python.md#rdkit) and [Gemmi](../interfaces/python.md#gemmi) objects,
  and get charges back as arrays without writing files.
- A [C++ library](../interfaces/cpp.md) with an installed CMake package.
- Charges for every conformer: all models of a PDB or mmCIF file and all conformers of an RDKit
  molecule are calculated, each with its own charges for methods that use geometry.
- [Explanations](../guide/methods.md#choosing-a-method) of why a method or parameter set does not apply,
  and a [result JSON](../guide/results.md) that records the model and settings of every calculation.
- [Cutoff and cover](../guide/large-structures.md) for ABEEM and the SQE methods, in addition to the
  EEM-type methods.
- [Fixed ions](../guide/fixed-ions.md) for structures with monatomic ions.
- Method implementations checked against their publications; several were corrected, so their charges
  can differ from ChargeFW2. The [method descriptions](../guide/methods.md#method-descriptions) give the
  equations ChargeFW uses.

## Moving from ChargeFW2

On the command line:

| ChargeFW2 | ChargeFW |
| --- | --- |
| `chargefw2 --mode info --input-file molecules.sdf` | `chargefw inspect molecules.sdf` |
| `chargefw2 --mode suitable-methods --input-file molecules.sdf` | `chargefw applicability molecules.sdf` |
| `chargefw2 --mode charges --input-file molecules.sdf --method eem --chg-out-dir out` | `chargefw calculate --method eem molecules.sdf out` |
| `--par-file EEM_00_NEEMP_ccd2016_npa` | `--parameter-set EEM_NEEMP_ccd2016_npa` |

In Python:

| ChargeFW2 | ChargeFW |
| --- | --- |
| `chargefw2.Molecules("molecules.sdf")` | `chargefw.io.read("molecules.sdf", format="sdf")` |
| `chargefw2.get_available_methods()` | `chargefw.methods` |
| `chargefw2.get_available_parameters("eem")` | `chargefw.methods["eem"].parameter_sets` |
| `chargefw2.get_suitable_methods(molecules)` | `chargefw.assess(molecules).plans` |
| `chargefw2.get_best_parameters(molecules, "eem")` | `chargefw.assess(molecules, method="eem").default_plan.parameter_set` |
| `chargefw2.calculate_charges(molecules, "eem", "EEM_00_NEEMP_ccd2016_npa")` | `chargefw.calculate(molecules, method="eem", parameter_set="EEM_NEEMP_ccd2016_npa")` |
| `chargefw2.save_charges(...)` | `chargefw.io.write(path, result, format=...)` |

Things to know when you move over:

- **Parameter-set IDs** no longer contain a number, so `EEM_00_NEEMP_ccd2016_npa` is now
  `EEM_NEEMP_ccd2016_npa`. The [parameter-set reference](../reference/parameter-sets.md) lists them.
- **Results** are assignments per molecule and conformer, in input order, not a dictionary keyed by
  molecule name; [Results](../guide/results.md) describes them.
- **Automatic selection** tries methods and parameter sets in a fixed priority order, the same for
  every structure. Name the method and parameter set to get a particular model.
- **One model per calculation:** all molecules of a call are calculated with the same model, and the
  calculation fails if it does not apply to one of them. Calculate molecules
  [one by one](../interfaces/python.md#one-model-for-all-molecules) to let each get its own model.
- **Structure selection** replaces the HETATM and water options: `polymers-and-ligands` leaves out
  water and `polymers` leaves out ligands too; see [Reading the
  input](../interfaces/cli.md#reading-the-input).
- **Bonds of ligands:** residue templates cover common amino acids, nucleotides, and water. Ligands in
  PDB files need `CONECT` records; mmCIF files downloaded from the PDB contain the bonds of their
  ligands.
- **Output** is result JSON, MOL2, or mmCIF. For a table of charges, use the
  [CSV example](../examples/calculate_file.md).
- **Parameter files of your own** can be loaded from C++; the command line and Python use the installed
  parameter sets.
- **SMP/QEq** has no cutoff or cover; calculate it in full.
