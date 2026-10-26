# ChargeFW for Python

ChargeFW calculates empirical partial atomic charges for small molecules and biomolecular structures.

ChargeFW is the new implementation of ChargeFW2, the calculation core of
[Atomic Charge Calculator III](https://acc.biodata.ceitec.cz/)
([Raček et al., 2026](https://doi.org/10.1093/nar/gkag379)). It calculates charges with the 20 methods
and the parameter sets of ChargeFW2, works with molecules in memory and with RDKit and Gemmi objects,
gives charges for every conformer, explains why a method does not apply, and records every calculation.
[Differences from ChargeFW2](https://chargefw.readthedocs.io/en/latest/about/chargefw2.html) describes
what changed and how to move over.

## Installation

ChargeFW 3 is currently a pre-release, so pip needs `--pre`. It requires Python 3.11 or newer on Linux
x86-64:

```bash
pip install --pre chargefw
```

The RDKit and Gemmi integrations are optional:

```bash
pip install --pre "chargefw[rdkit]"
pip install --pre "chargefw[gemmi]"
```

## Usage

```python
import chargefw

molecules = chargefw.io.read("molecule.sdf", format="sdf")
result = chargefw.calculate(molecules, method="eem", parameter_set="EEM_Cheminf_b3lyp_npa")
print(result.assignments[0].values)

chargefw.io.write("charges.json", result, format="result-json")
```

Molecules can also be built from NumPy arrays:

```python
molecule = chargefw.Molecule(
    atomic_numbers=[8, 1, 1],
    bonds=[[0, 1, 1], [0, 2, 1]],
    coordinates=[[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]],
)
result = chargefw.calculate(molecule, method="eem")
```

Name both the method and the parameter set when the results must be reproducible.

## Documentation

The documentation is at [chargefw.readthedocs.io](https://chargefw.readthedocs.io/):

- [First calculation](https://chargefw.readthedocs.io/en/latest/getting-started/first-calculation.html):
  an example with its output.
- [Preparing input](https://chargefw.readthedocs.io/en/latest/guide/input.html): what a molecule needs
  before a calculation.
- [Methods and parameter sets](https://chargefw.readthedocs.io/en/latest/guide/methods.html): the
  available methods and how to choose one.
- [Python](https://chargefw.readthedocs.io/en/latest/interfaces/python.html): the complete Python
  interface.
- [Examples](https://chargefw.readthedocs.io/en/latest/examples/index.html): complete scripts.

The source code is at [github.com/sb-ncbr/chargefw](https://github.com/sb-ncbr/chargefw). ChargeFW is
distributed under the MIT license.
