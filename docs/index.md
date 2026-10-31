# ChargeFW

ChargeFW calculates empirical partial atomic charges for small molecules and biomolecular structures.
It provides a Python package, a command-line tool, and a C++ library, all using the same calculation
engine, methods, and bundled parameter sets.

ChargeFW is the new implementation of ChargeFW2, the calculation core of
[Atomic Charge Calculator III](https://acc.biodata.ceitec.cz/)
([Raček et al., 2026](https://doi.org/10.1093/nar/gkag379)). It has the same 20 methods and parameter
sets, and adds:

- a [Python package](interfaces/python.md) that works with molecules in memory, including RDKit and
  Gemmi objects;
- [charges for every conformer](guide/results.md) of a molecule and every model of a structure;
- [reasons](guide/methods.md#choosing-a-method) why a method or parameter set does not apply;
- a [record](guide/results.md#what-the-result-records) of the method, parameter set, and settings of
  every calculation;
- [faster approximate calculation](guide/large-structures.md) of large structures for more methods;
- [fixed ions](guide/fixed-ions.md) for structures with metal and other monatomic ions.

[Differences from ChargeFW2](about/chargefw2.md) describes what changed and how to move over.

You give ChargeFW a molecule with its bonds, formal charges, and coordinates where the method needs
them. It does not change your structure; [Preparing input](guide/input.md) explains what to check
before a calculation.

## Installation

```bash
pip install --pre chargefw
```

The [installation page](getting-started/installation.md) also covers the Docker image and building
from source. [First calculation](getting-started/first-calculation.md) then takes you through one
example on the command line and in Python.

```{toctree}
:hidden:
:caption: Getting started

getting-started/installation
getting-started/first-calculation
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
reference/python-api
```

```{toctree}
:hidden:
:caption: About

about/versions
about/chargefw2
about/citing
```

ChargeFW is developed at [sb-ncbr/chargefw](https://github.com/sb-ncbr/chargefw) and distributed
under the MIT license.
