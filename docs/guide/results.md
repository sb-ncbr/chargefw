# Results

A calculation returns the charges of every atom, in the order of the atoms in your input, together with
a record of the model that produced them.

## Charges per molecule and conformer

Methods that use 3D coordinates give each conformer of a molecule its own charges. The models of a PDB
or mmCIF file, the conformers of an RDKit molecule, and the conformers in molecule JSON all count as
conformers. Methods that use only the molecular graph give one set of charges per molecule.

Each set of charges is an assignment. Here ethanol has three RDKit conformers:

```python
for method in ("eem", "peoe"):
    result = chargefw.calculate(molecule, method=method)
    for assignment in result.assignments:
        print(method, assignment.molecule_index, assignment.conformer_index,
              assignment.values[:3].round(3))
```

```text
eem 0 0 [-0.503 -0.017 -0.876]
eem 0 1 [-0.488 -0.042 -0.867]
eem 0 2 [-0.488 -0.021 -0.876]
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

`assignment.atom_ids` holds one ID per atom: its ID in the input file, or its position in the file when
the format has no atom IDs. For PDB and mmCIF files with several models, these IDs are those of the
first model. The atom sites of the model that a conformer came from are in
`molecule.source_mapping.conformers[assignment.conformer_index].sites`, each with its `id`. mmCIF
output and attaching charges to a Gemmi document use this mapping for you.

## What the result records

Every successful result records the model that was used: the method, the parameter set, the method
options with their values, and the execution mode with its radius. In Python, `result.plan` holds
them; in result JSON, `calculation_provenance.effective` does. Use these values to repeat a calculation
exactly, or to report it.

Result JSON also records what you requested, any warnings, the time each step took, and the peak memory.
The [result JSON reference](../reference/formats.md#chargefw-result-json-10) lists every field.

## When a calculation fails

A calculation fails without charges when no method and parameter set apply to your molecules, when the
calculation does not converge, or when it is cancelled. The command line still writes the result JSON,
with the reasons, and exits with a [status](../interfaces/cli.md#exit-statuses) other than 0. Python
raises a `chargefw.ChargeFWError`; its `result` attribute holds the same information. Here the
parameter set has no values for chlorine, and `molecules` holds chlorobenzene:

```python
try:
    result = chargefw.calculate(molecules, method="eem", parameter_set="EEM_Baek1991")
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
| Python objects | Further calculation: charges as NumPy arrays, attached to RDKit molecules or Gemmi documents |

MOL2 and mmCIF files contain charges rounded to four decimal places. The command line writes them with
`--output-mol2` and `--output-mmcif`; in Python, use `chargefw.io.write()`. The
[format reference](../reference/formats.md#charge-output) describes each format.
