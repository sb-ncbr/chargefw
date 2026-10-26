# Versions and changes

The native library, the command-line application, and the Python package share one
[PEP 440](https://peps.python.org/pep-0440/) version, which continues ChargeFW2's lineage. Pre-releases
such as `3.0.0a1` may change any interface. From 3.0.0, incompatible changes to documented interfaces
wait for a new major version, and every change that affects calculated charges is listed in the
[changelog](https://github.com/sb-ncbr/chargefw/blob/main/CHANGELOG.md). New parameter sets can change
automatic selection, so name the method and parameter set explicitly when results must stay the same
across releases.
