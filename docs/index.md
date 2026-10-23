# ChargeFW

ChargeFW calculates empirical partial atomic charges. It provides a toolkit-neutral C++23 library,
a command-line application, and Python bindings over the same calculation engine, with bundled
methods and parameter sets, deterministic method selection, and source-ordered results.

ChargeFW expects a molecular graph, formal charges, and any coordinates required by the selected
method; it does not prepare structures. [Project design](PROJECT.md) describes the scientific scope.

## Installation

- **Python:** `pip install chargefw` (Python 3.11 or newer); see [Python package](PYTHON.md).
- **Command line and C++:** build from source as described in [Native C++ library](NATIVE.md).

```{toctree}
:maxdepth: 2
:caption: User guide

PYTHON
recipes/README
CLI
FORMATS
PARAMETERS
NATIVE
PROJECT
```

ChargeFW is developed at [sb-ncbr/chargefw](https://github.com/sb-ncbr/chargefw) and distributed
under the MIT license.
