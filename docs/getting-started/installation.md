# Installation

## Python package

You need Python 3.12 or newer on Linux (x86-64 or aarch64), macOS 14 or newer on Apple silicon, or
64-bit Windows. ChargeFW 3 is currently a pre-release, so tell pip to accept it with `--pre`:

```bash
pip install --pre chargefw
```

Without `--pre`, pip installs 0.1.7, the last release of the previous series. One wheel serves CPython
3.12 and newer, and pip installs NumPy 2.0 or newer with it.

On Linux, the wheels need glibc 2.27 or newer (2.26 on aarch64), as in Ubuntu 18.04, Debian 10, RHEL 8,
and later releases. There are no wheels for Alpine Linux and other systems based on musl; pip then tries
to build ChargeFW from source, which is not supported. Use the [Docker image](#command-line-in-docker)
for the command line there.

On Python 3.11, pip reports no error but installs an older ChargeFW release. If your system Python is
3.11, as on Debian 12, create an environment with a newer Python. [uv](https://docs.astral.sh/uv/)
downloads one for you:

```bash
uv venv --python 3.12
source .venv/bin/activate
uv pip install --pre chargefw
```

Reading and writing Gemmi structures, and converting RDKit molecules, need the optional extras:

```bash
pip install --pre "chargefw[gemmi]"
pip install --pre "chargefw[rdkit]"
```

The `gemmi` extra installs exactly Gemmi 0.7.4 and replaces any other Gemmi version in the
environment.

pip also installs the `chargefw` command line. Check both:

```bash
python -c "import chargefw; print(chargefw.__version__)"
chargefw --help
```

## Command line in Docker

The `sbncbr/chargefw` image on Docker Hub contains the `chargefw` command and the bundled parameter
sets. Mount the directory with your structures and run a command in it:

```bash
docker run --rm --user "$(id -u):$(id -g)" -v "$PWD:/work" -w /work \
    sbncbr/chargefw:3.0.0rc2 calculate molecule.sdf output
```

The charges are written to `output/`. `--user` makes the output files yours rather than root's. Every
release has its own tag; from 3.0.0 on, `latest` points to the newest final release.

## From source

The C++ library and the source build of the command line are supported on Linux. To use the library,
build and install ChargeFW as described in
[Build and link](../interfaces/cpp.md#build-and-link). Set `CHARGEFW_BUILD_CLI=ON` to build the command
line as well. The installed `chargefw` finds its parameter sets relative to the installed library, so it
works from any directory:

```bash
_install/bin/chargefw --help
```
