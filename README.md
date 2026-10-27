# ChargeFW

ChargeFW calculates empirical partial atomic charges for small molecules and biomolecular structures,
from Python, from the command line, or from C++.

ChargeFW is the new implementation of ChargeFW2, the calculation core of
[Atomic Charge Calculator III](https://acc.biodata.ceitec.cz/)
([Raček et al., 2026](https://doi.org/10.1093/nar/gkag379)). It calculates charges with the 20 methods
and the parameter sets of ChargeFW2, and adds an in-memory Python interface with RDKit and Gemmi
integration, a C++ library, charges for every conformer, explanations of why a method does not apply,
a record of every calculation, cutoff and cover for more methods, and fixed ions. Method implementations
were checked against their publications. [Differences from
ChargeFW2](https://chargefw.readthedocs.io/en/latest/about/chargefw2.html) describes what changed and how
to move over.

## Installation

```bash
pip install --pre chargefw
```

ChargeFW 3 is currently a pre-release, so pip needs `--pre`. pip installs the Python package and the
`chargefw` command line. The command line is also available as the Docker image `sbncbr/chargefw`, and
everything can be built from source; see
[Installation](https://chargefw.readthedocs.io/en/latest/getting-started/installation.html).

## Quick start

In Python:

```python
import chargefw

molecules = chargefw.io.read("molecule.sdf", format="sdf")
result = chargefw.calculate(molecules, method="eem", parameter_set="EEM_Cheminf_b3lyp_npa")
print(result.assignments[0].values)
```

On the command line:

```bash
chargefw calculate --method eem --parameter-set EEM_Cheminf_b3lyp_npa molecule.sdf output
```

or with the Docker image:

```bash
docker run --rm --user "$(id -u):$(id -g)" -v "$PWD:/work" -w /work sbncbr/chargefw:3.0.0a4 \
    calculate --method eem --parameter-set EEM_Cheminf_b3lyp_npa molecule.sdf output
```

[First calculation](https://chargefw.readthedocs.io/en/latest/getting-started/first-calculation.html)
walks through an example with its output.

## Documentation

The documentation is at [chargefw.readthedocs.io](https://chargefw.readthedocs.io/):

- [Preparing input](https://chargefw.readthedocs.io/en/latest/guide/input.html): what a molecule needs
  before a calculation.
- [Methods and parameter sets](https://chargefw.readthedocs.io/en/latest/guide/methods.html): the
  available methods, their equations, and how to choose one.
- [Python](https://chargefw.readthedocs.io/en/latest/interfaces/python.html),
  [command line](https://chargefw.readthedocs.io/en/latest/interfaces/cli.html), and
  [C++](https://chargefw.readthedocs.io/en/latest/interfaces/cpp.html).
- [Examples](https://chargefw.readthedocs.io/en/latest/examples/index.html): complete scripts for files,
  RDKit, Gemmi, and analysis.

## Development

- [Changelog](CHANGELOG.md): changes in each release.
- [Development guide](DEVELOPMENT.md): building, testing, and validating the repository.
- [Release guide](RELEASING.md): publishing a release.
- [Unfinished work](TODO.md): current gaps.
- [Repository change policy](AGENTS.md): architecture and validation rules.

To cite ChargeFW, use the metadata in [CITATION.cff](CITATION.cff). ChargeFW is distributed under the
MIT license.
