# Installation

## Python package

You need Python 3.11 or newer on Linux x86-64. ChargeFW 3 is currently a pre-release, so tell pip to
accept it with `--pre`:

```bash
pip install --pre chargefw
```

Without `--pre`, pip installs 0.1.7, the last release of the previous series. Wheels are published for
CPython 3.11 to 3.14, and pip installs NumPy 2.0 or newer with them.

Reading and writing Gemmi structures, and converting RDKit molecules, need the optional extras:

```bash
pip install --pre "chargefw[gemmi]"
pip install --pre "chargefw[rdkit]"
```

Check the installation:

```bash
python -c "import chargefw; print(chargefw.__version__)"
```

## Command line in Docker

The `sbncbr/chargefw` image on Docker Hub contains the `chargefw` command and the bundled parameter
sets. Mount the directory with your structures and run a command in it:

```bash
docker run --rm --user "$(id -u):$(id -g)" -v "$PWD:/work" -w /work \
    sbncbr/chargefw:3.0.0a3 calculate molecule.sdf output
```

The charges are written to `output/`. `--user` makes the output files yours rather than root's. Every
release has its own tag; from 3.0.0 on, `latest` points to the newest final release.

## From source

To use the C++ library, or the command line without Docker, build and install ChargeFW as described in
[Build and link](../interfaces/cpp.md#build-and-link). Set `CHARGEFW_BUILD_CLI=ON` to build the command
line as well. The installed `chargefw` finds its parameter sets relative to the installed library, so it
works from any directory:

```bash
_install/bin/chargefw --help
```
