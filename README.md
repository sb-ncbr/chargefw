# ChargeFW

[![CI](https://github.com/sb-ncbr/chargefw/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/sb-ncbr/chargefw/actions/workflows/ci.yml)
[![Documentation](https://img.shields.io/readthedocs/chargefw/stable?label=docs)](https://chargefw.readthedocs.io/en/stable/)
[![PyPI](https://img.shields.io/pypi/v/chargefw)](https://pypi.org/project/chargefw/)
[![Cite: NAR 2026](https://img.shields.io/badge/cite-NAR%202026-blue)](https://doi.org/10.1093/nar/gkag379)

ChargeFW calculates empirical partial atomic charges for small molecules and biomolecular structures,
from Python, from the command line, or from C++.

ChargeFW is the new implementation of ChargeFW2, the calculation core of
[Atomic Charge Calculator III](https://acc.biodata.ceitec.cz/)
([Raček et al., 2026](https://doi.org/10.1093/nar/gkag379)). It has the same 20 methods, each
implementation checked against its publication, and the same parameter sets, and adds:

- a Python package that works with molecules in memory, including RDKit and Gemmi objects;
- charges for every conformer of a molecule and every model of a structure;
- reasons why a method or parameter set does not apply;
- a record of the method, parameter set, and settings of every calculation;
- faster approximate calculation of large structures for more methods;
- fixed ions for structures with metal and other monatomic ions.

[Differences from ChargeFW2](https://chargefw.readthedocs.io/en/stable/about/chargefw2.html)
describes what changed and how to move over.

## Installation

```bash
pip install chargefw
```

pip installs the Python package and the `chargefw` command line. The command line is also available as
the Docker image `sbncbr/chargefw`, and everything can be built from source; see
[Installation](https://chargefw.readthedocs.io/en/stable/getting-started/installation.html).

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
docker run --rm --user "$(id -u):$(id -g)" -v "$PWD:/work" -w /work sbncbr/chargefw \
    calculate --method eem --parameter-set EEM_Cheminf_b3lyp_npa molecule.sdf output
```

[First calculation](https://chargefw.readthedocs.io/en/stable/getting-started/first-calculation.html)
walks through an example with its output.

## Documentation

The documentation is at [chargefw.readthedocs.io](https://chargefw.readthedocs.io/):

- [Preparing input](https://chargefw.readthedocs.io/en/stable/guide/input.html): what a molecule needs
  before a calculation.
- [Methods and parameter sets](https://chargefw.readthedocs.io/en/stable/guide/methods.html): the
  available methods, their equations, and how to choose one.
- [Python](https://chargefw.readthedocs.io/en/stable/interfaces/python.html),
  [command line](https://chargefw.readthedocs.io/en/stable/interfaces/cli.html), and
  [C++](https://chargefw.readthedocs.io/en/stable/interfaces/cpp.html).
- [Examples](https://chargefw.readthedocs.io/en/stable/examples/index.html): complete scripts for files,
  RDKit, Gemmi, and analysis.

## Development

- [Changelog](CHANGELOG.md): changes in each release.
- [Contributing](CONTRIBUTING.md): reporting a problem and preparing a pull request.
- [Development guide](DEVELOPMENT.md): building, testing, and validating the repository.
- [Release guide](RELEASING.md): publishing a release.
- [Unfinished work](TODO.md): current gaps.
- [Repository change policy](AGENTS.md): architecture and validation rules.

To cite ChargeFW, use the metadata in [CITATION.cff](CITATION.cff). ChargeFW is distributed under the
MIT license.
