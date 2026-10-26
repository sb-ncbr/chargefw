# ChargeFW for Python

ChargeFW calculates empirical partial atomic charges through Python bindings to its native C++
calculation engine. It provides bundled methods and parameter sets, deterministic applicability and
execution planning, source-ordered result mapping, molecular-file input, result JSON, MOL2 and mmCIF
output, and optional Gemmi and RDKit integration.

ChargeFW expects prepared molecular graphs with formal charges and any coordinates required by the
selected method. It does not parse SMILES, add hydrogens, choose protonation states, perceive arbitrary
bonds, or generate coordinates.

## Installation

ChargeFW requires Python 3.11 or newer:

```bash
pip install chargefw
```

Install an optional toolkit integration when converting existing toolkit objects:

```bash
pip install "chargefw[gemmi]"
pip install "chargefw[rdkit]"
```

The base package can read PDB and mmCIF text through its compiled Gemmi dependency without installing
the Gemmi Python package.

## Python usage

### From a molecular file

```python
import chargefw

molecules = chargefw.io.read("molecule.sdf", format="sdf")
result = chargefw.calculate(
    molecules,
    method="qeq",
    parameter_set="QEq_original",
    execution="full",
)
print(result.assignments[0].values)
```

For reproducible scientific work, select both the method and its parameter set when the method uses
one. Automatic selection is deterministic, but catalog priority is not a recommendation that a model
or parameterization is scientifically preferable for a particular molecule.

### From arrays

```python
import chargefw

molecule = chargefw.Molecule(
    atomic_numbers=[8, 1, 1],
    formal_charges=[0, 0, 0],
    bonds=[[0, 1, 1], [0, 2, 1]],
    coordinates=[
        [0.0, 0.0, 0.0],
        [0.96, 0.0, 0.0],
        [-0.24, 0.93, 0.0],
    ],
    name="water",
)

result = chargefw.calculate(
    molecule,
    method="qeq",
    parameter_set="QEq_original",
    execution="full",
)
print(result.assignments[0].values)
```

Write a complete result record with `chargefw.io.write()`:

```python
chargefw.io.write("charges.json", result, format="result-json")
```

## Documentation

The full documentation is at [chargefw.readthedocs.io](https://chargefw.readthedocs.io/):

- **[Python API](https://chargefw.readthedocs.io/en/latest/interfaces/python.html)** - read molecules,
  calculate charges, and inspect methods and results.
- **[Examples](https://chargefw.readthedocs.io/en/latest/examples/index.html)** - run complete
  SDF, Gemmi, RDKit, and parameter-set workflows.
- **[Molecular formats](https://chargefw.readthedocs.io/en/latest/reference/formats.html)** - check supported
  input, output, parsing, and mapping semantics.
- **[Methods and parameter sets](https://chargefw.readthedocs.io/en/latest/guide/methods.html)** - see
  the available methods and how they are selected.
- **[Source repository](https://github.com/sb-ncbr/chargefw)** - browse the source and report
  issues.

ChargeFW is distributed under the MIT license.
