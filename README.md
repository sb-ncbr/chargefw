# ChargeFW

ChargeFW is a C++23 framework for empirical partial atomic-charge calculation. It provides a
toolkit-neutral native library, a command-line application, and Python bindings over the same
assessment and calculation engine.

The current implementation includes bundled methods and parameter sets, deterministic applicability and
execution planning, source-ordered result mapping, molecular-file adapters, and full or explicit reduced
calculations.

ChargeFW expects a molecular graph, formal charges, and any coordinates required by the selected
method. It does not parse SMILES, add hydrogens, assign protonation states, perceive arbitrary bonds,
or generate coordinates.

## Quick start

Build the [native library](docs/NATIVE.md#build-and-link) with
`CHARGEFW_BUILD_CLI=ON`, then calculate charges for an SDF file:

```bash
_install/bin/chargefw calculate molecule.sdf output
```

For Python installation and package status, see the [Python package guide](docs/PYTHON.md#installation).

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

## Documentation

The documentation is published at [chargefw.readthedocs.io](https://chargefw.readthedocs.io/). Its
sources are in this repository:

- **[Python API](docs/PYTHON.md)** - read molecules, calculate charges, and inspect methods and results.
- **[Python recipes](docs/recipes/README.md)** - run complete SDF, Gemmi, RDKit, and parameter-set
  workflows.
- **[Command-line interface](docs/CLI.md)** - calculate and inspect molecular files from the CLI.
- **[Molecular formats](docs/FORMATS.md)** - check supported input, output, parsing, and mapping
  semantics.
- **[Parameter sets](docs/PARAMETERS.md)** - understand parameter data and matching behavior.
- **[Native C++ API](docs/NATIVE.md)** - embed ChargeFW in a C++ application.
- **[Project design](docs/PROJECT.md)** - understand the architecture and scientific scope.

## Development

- **[Changelog](CHANGELOG.md)** - see changes in each release.
- **[Development guide](DEVELOPMENT.md)** - configure, build, test, and validate the repository.
- **[Release guide](RELEASING.md)** - build, verify, and publish Python distributions.
- **[Unfinished work](TODO.md)** - see current project gaps.
- **[Repository change policy](AGENTS.md)** - preserve architecture and choose proportionate validation.

ChargeFW is distributed under the MIT license.
