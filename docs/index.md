# ChargeFW

ChargeFW calculates empirical partial atomic charges for small molecules and biomolecular structures.
It provides a Python package, a command-line tool, and a C++ library, all using the same calculation
engine, methods, and bundled parameter sets.

You give ChargeFW a molecule with its bonds, formal charges, and coordinates where the method needs
them. It does not change your structure; [Preparing input](guide/input.md) explains what to check
before a calculation.

## Installation

```bash
pip install --pre chargefw
```

The [installation page](getting-started/installation.md) also covers the Docker image and building
from source.

```{toctree}
:hidden:
:caption: Getting started

getting-started/installation
```

```{toctree}
:hidden:
:caption: User guide

guide/input
guide/methods
guide/large-structures
guide/fixed-ions
guide/results
```

```{toctree}
:hidden:
:caption: Interfaces

interfaces/python
interfaces/cli
interfaces/cpp
```

```{toctree}
:hidden:
:caption: Examples

examples/index
```

```{toctree}
:hidden:
:caption: Reference

reference/formats
reference/parameter-sets
reference/parameter-format
```

```{toctree}
:hidden:
:caption: About

about/versions
```

ChargeFW is developed at [sb-ncbr/chargefw](https://github.com/sb-ncbr/chargefw) and distributed
under the MIT license.
