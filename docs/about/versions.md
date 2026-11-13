# Versions and changes

The library, the command line, and the Python package share one version, such as `3.0.0`. A patch
release, such as `3.0.1`, only fixes bugs. A minor release, such as `3.1.0`, may also change
documented interfaces; the [changelog](https://github.com/sb-ncbr/chargefw/blob/main/CHANGELOG.md)
marks every such change and says what to do instead. To stay within one minor version, install
`chargefw~=3.0.0`. `chargefw.parameterization` is experimental and may change in any release.
Rebuild C++ applications against every release.

A bundled parameter set never changes its values under the same ID, and every change to calculated
charges is listed in the changelog. Automatic selection may choose differently in a minor release,
so name the method, the parameter set, and the execution when results must stay the same.
