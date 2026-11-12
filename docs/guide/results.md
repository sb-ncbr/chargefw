# Results

A calculation returns the charges of every atom, in the order of the atoms in your input, together with
a record of the model that produced them.

## Charges per molecule and conformer

Methods that use 3D coordinates give each conformer of a molecule its own charges. The models of a PDB
or mmCIF file, the conformers of an RDKit molecule, and the conformers in molecule JSON all count as
conformers. Methods that use only the molecular graph give one set of charges per molecule.

ChargeFW calls each set of charges an assignment. Here ethanol gets three RDKit conformers; the
RDKit integration needs `pip install --pre "chargefw[rdkit]"`:

```python
from rdkit import Chem
from rdkit.Chem import AllChem

import chargefw
import chargefw.io.rdkit

ethanol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
AllChem.EmbedMultipleConfs(ethanol, numConfs=3, randomSeed=42)
molecule = chargefw.io.rdkit.from_mol(ethanol)

for method in ("eem", "peoe"):
    result = chargefw.calculate(molecule, method=method)
    for assignment in result.assignments:
        print(
            method,
            assignment.molecule_index,
            assignment.conformer_index,
            assignment.values[:3].round(3),
        )
```

```text
eem 0 0 [-0.476 -0.015 -0.878]
eem 0 1 [-0.492 -0.028 -0.868]
eem 0 2 [-0.486 -0.01  -0.875]
peoe 0 None [-0.042  0.041 -0.395]
```

EEM depends on geometry, so the charges differ between conformers. PEOE does not, and its single
assignment has no conformer index. To get one assignment directly, use
`result.assignment(molecule=0, conformer=1)`; omit `conformer` for methods that do not use coordinates.

The result JSON written by the command line lists the same assignments under `results[i].assignments`,
each with a `target` naming its molecule and conformer.

## Atom order

The charges of an assignment follow the atoms of the molecule as ChargeFW read it: the order of the
input file, or of the RDKit molecule or arrays you passed in. For PDB and mmCIF input, only the atoms
that pass the [structural selection](../reference/formats.md#structural-selection) are included, still
in file order.

`assignment.atom_ids` holds one ID per atom: its ID in the input file, or its position in the file
when the format has no atom IDs. For PDB and mmCIF files with several models, these IDs are those of
the first model; the atoms of the other models are listed, with their IDs, in
`molecule.source_mapping.conformers[assignment.conformer_index].sites`. mmCIF output and
`attach_charges()` for Gemmi documents put each charge on the right atom of every model for you.

## What the result records

Every successful result records the model that was used: the method, the parameter set, the method
options with their values, and the execution mode with its radius. In Python, `result.plan` holds
them; in result JSON, `calculation_provenance.effective` does. Use these values to repeat a calculation
exactly, or to report it.

Result JSON also records what you requested, any warnings, the time each step took, and the peak memory.
The [result JSON reference](../reference/formats.md#chargefw-result-json-10) lists every field.

## Warnings

ChargeFW calculates charges also when it finds a problem with the input or with the charges, and
reports the problem as a warning. Read the warnings: some of them mean that the charges are probably
wrong.

On the command line, each warning is printed on standard error with the molecule it concerns;
warnings alone do not change the exit status:

```console
$ chargefw calculate --method eem 1CRN.cif output
Warning: molecule 1 ('1CRN'): the molecule contains carbon but no hydrogen atoms; if its hydrogens were not added, add them before calculating
Wrote output/1CRN.chargefw.json
```

When several molecules of a file get the same warning, it is printed once, with their number and the
first of them, such as `Warning: 2 of 3 molecules: hydrogens_missing, for example molecule 1 …`; the
result JSON lists each one.

Result JSON records each warning with its code in the `diagnostics` of the molecule's record, or of
the document when it concerns no single molecule. In Python, the warnings that can make the charges
wrong are also [Python warnings](../interfaces/python.md#warnings), which Python prints when your
script runs; every warning is recorded in `molecule.import_issues`, `result.input_warnings`,
`result.charge_warnings`, or `result.plan.warnings`.

| Code | Meaning | What to do | Python warning |
| --- | --- | --- | --- |
| `charge_out_of_range` | Charges are more than 6 e away from the formal charges; they are probably wrong ([below](#charges-that-are-probably-wrong)) | Check the input, or use another method or parameter set | `ChargeWarning` |
| `atoms_too_close` | Atoms of a conformer overlap, which makes the charges of methods that use coordinates wrong ([Preparing input](input.md)) | Check the coordinates of the named atoms | `InputWarning` |
| `hydrogens_missing` | The molecule contains carbon but no hydrogen atoms | Check whether hydrogens are missing; add them if needed | `InputWarning` |
| `unbonded_residue_atoms` | Atoms of a residue in a PDB or mmCIF file have no bonds | Add their bonds to the file | `InputWarning` |
| `symmetry_connections_ignored` | Bonds to symmetry copies of atoms in a crystal structure were left out | If the copies belong to your molecule, build them before calculating | `InputWarning` |
| `aromatic_bonds_converted` | Aromatic bonds were read as single bonds, which changes the charges of methods that use bond orders | Write the file with Kekulé bonds |  |
| `coordinates_2d_ignored` | A 2D drawing was read without coordinates | For methods that use coordinates, generate 3D coordinates |  |
| `partial_charges_ignored` | The file contains partial charges, which ChargeFW does not use | Nothing |  |
| `v3000_bond_order_ignored` | Coordinate or hydrogen bonds of a V3000 record were left out | Usually nothing |  |
| `resource_threshold_exceeded` | You requested full execution of a structure above the [size limit](large-structures.md) | Nothing; expect a long calculation |  |

## Charges that are probably wrong

ChargeFW compares every calculated charge with the formal charge of its atom. When a charge is more
than 6 e away from it, the molecule gets a `charge_out_of_range` warning, or the conformer does for
methods that use coordinates. Recommended methods give charges within about 5 e of the formal charges
for correctly prepared molecules, the phosphorus of phosphates included. Larger charges come from
atoms that overlap, from structures with errors, or from a method that does not work for the
molecule, such as {ref}`GDAC <gdac>` for amides and peptides. ChargeFW still returns the
charges. Check the input,
as described in [Preparing input](input.md), or use another method or parameter set.

The warning names the atom furthest from its formal charge and says how many atoms are affected. In
Python, `result.charge_warnings` lists the warnings, and the `atom_indices` of each has every
affected atom. In result JSON, the warning in the molecule's record has the same `atom_indices`, so
that a viewer can select those atoms. [Fixed ions](fixed-ions.md) are not checked, because their
charge is set.

## When a calculation fails

A calculation fails without charges when no method and parameter set apply to your molecules, when the
calculation does not converge, or when it is cancelled. With SQE, SQE+q0, and SQE+qp, it also fails
when the equations for the split charges cannot be factorized as positive definite. This happens when
the parameter set gives the molecule no single energy minimum at its geometry, and can also happen
through rounding when it comes close to that. With cutoff and cover execution, each fragment's
equations are checked, not those of the whole molecule. Check the geometry, or use another parameter
set. The command line still writes the result JSON, with the reasons, and exits with a
[status](../interfaces/cli.md#exit-statuses) other than 0. Python raises a `chargefw.ChargeFWError`;
its `result` attribute holds the same information. Here the parameter set has no values for chlorine
in chlorobenzene:

```python
chlorobenzene = Chem.AddHs(Chem.MolFromSmiles("Clc1ccccc1"))
AllChem.EmbedMolecule(chlorobenzene, randomSeed=42)
Chem.Kekulize(chlorobenzene, clearAromaticFlags=True)
molecule = chargefw.io.rdkit.from_mol(chlorobenzene)

try:
    result = chargefw.calculate(molecule, method="eem", parameter_set="EEM_Baek1991")
except chargefw.ChargeFWError as error:
    for rejection in error.result.rejections:
        for issue in rejection.issues:
            print(rejection.method.id, issue.message)
```

```text
eem molecule 1: parameter set 'EEM_Baek1991' has no atom parameter matching atom 1 (source name 'Cl1', Cl, formal charge 0)
```

[Choosing a method](methods.md#choosing-a-method) shows how to see the reasons before calculating.

## Choosing an output format

| Format | Use it for |
| --- | --- |
| Result JSON | Scripts and analysis: full-precision charges and the complete record of the calculation |
| MOL2 | Small molecules in molecular viewers and docking programs |
| mmCIF | Proteins and other structures: charges stay linked to chains, residues, and atoms |
| PQR | Electrostatics of proteins read from PDB or mmCIF, for example with APBS; one model only |
| Python objects | Further calculation: charges as NumPy arrays, attached to RDKit molecules or Gemmi documents |

MOL2, mmCIF, and PQR files contain charges rounded to four decimal places, and need coordinates. The
command line writes them with `--output-mol2`, `--output-mmcif`, and `--output-pqr`; in Python, use
`chargefw.io.write()`. PQR holds one model, so read a file with several models with
`--conformers first` on the command line, or `conformers="first"` in Python. The
[format reference](../reference/formats.md#charge-output) describes each format.
