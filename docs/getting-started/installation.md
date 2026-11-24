# Installation

## Python package

You need Python 3.12 or newer on Linux (x86-64 or aarch64), macOS 14 or newer on Apple silicon, or
64-bit Windows:

```bash
pip install chargefw
```

One wheel serves CPython 3.12 and newer, and pip installs NumPy 2.0 or newer with it.

On Linux, the wheels need glibc 2.27 or newer (2.26 on aarch64), as in Ubuntu 18.04, Debian 10, RHEL 8,
and later releases. There are no wheels for Alpine Linux and other systems based on musl; pip then tries
to build ChargeFW from source, which is not supported. Use the [Docker image](#command-line-in-docker)
for the command line there.

On x86-64, the wheels and the Docker image need a CPU with AVX2 and FMA (the x86-64-v3 level): Intel
Haswell (2013), AMD Excavator (2015), or newer. On an older CPU, `import chargefw` raises an
`ImportError`, and the Docker image stops with an illegal-instruction error. Build the library and
command line [from source](#from-source) with `CHARGEFW_ENABLE_NATIVE_OPTIMIZATIONS=ON` instead.

On Python 3.11, pip reports no error but installs an older ChargeFW release. If your system Python is
3.11, as on Debian 12, create an environment with a newer Python. [uv](https://docs.astral.sh/uv/)
downloads one for you:

```bash
uv venv --python 3.12
source .venv/bin/activate
uv pip install chargefw
```

Reading and writing Gemmi structures, and converting RDKit molecules, need the optional extras:

```bash
pip install "chargefw[gemmi]"
pip install "chargefw[rdkit]"
```

The `gemmi` extra installs exactly Gemmi 0.7.4 and replaces any other Gemmi version in the
environment.

The `ccd` extra installs `chargefw-ccd`, the bonds of the wwPDB Chemical Component Dictionary, which
ChargeFW uses for
[ligands and modified residues](../guide/input.md#bonds-of-ligands-and-modified-residues) when you
name it:

```bash
pip install "chargefw[ccd]"
```

Gemmi 0.7.4 has no wheel for Python 3.15, so pip builds it from source there, and that build fails
with nanobind 3. On Python 3.15, keep nanobind below version 3 for the build. This needs pip 25.3 or
newer, a C++ compiler, and the Python development headers (`python3-devel` on Fedora,
`python3-dev` on Debian and Ubuntu), and the build takes a few minutes:

```bash
echo "nanobind<3" > build-constraints.txt
pip install --build-constraint build-constraints.txt "chargefw[gemmi]"
```

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
    sbncbr/chargefw calculate molecule.sdf output
```

The charges are written to `output/`. `--user` makes the output files yours rather than root's. Without
a tag, Docker uses `latest`, the newest final release. Every release also has its own tag, such as
`sbncbr/chargefw:3.0.0`; name one when the results must stay the same.

## From source

The C++ library and the source build of the command line are supported on Linux. To use the library,
build and install ChargeFW as described in
[Build and link](../interfaces/cpp.md#build-and-link). Set `CHARGEFW_BUILD_CLI=ON` to build the command
line as well. The installed `chargefw` finds its parameter sets relative to the installed library, so it
works from any directory:

```bash
_install/bin/chargefw --help
```
