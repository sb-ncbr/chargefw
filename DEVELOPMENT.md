# ChargeFW development

This is the practical guide to building and testing ChargeFW. Repository policy and architecture rules
are in [AGENTS.md](AGENTS.md); public behavior is documented under [`docs/`](docs/).

## Requirements

Native development requires CMake 3.28 or newer, Ninja, and a C++23-capable GCC or Clang toolchain.
The checked-in presets set `CHARGEFW_FETCH_DEPENDENCIES=ON`, so CMake downloads pinned,
checksum-verified dependency sources. CI and the compatibility containers build against distribution
packages with `CHARGEFW_FETCH_DEPENDENCIES=OFF`; pass it to a preset to use installed packages locally.
The supported versions are listed in [docs/NATIVE.md](docs/NATIVE.md#build-and-link).

The checked-in presets also build the Python bindings. Their Python interpreter needs Python 3.11 or
newer, development headers, NumPy, Nanobind, Gemmi 0.7.4, and the test tools pytest, jsonschema, and
mypy. RDKit is optional; its tests are skipped when it is absent. The presets use `/usr/bin/python3`;
select another interpreter explicitly when needed:

```bash
cmake --preset gcc-debug \
    -DCHARGEFW_PYTHON_EXECUTABLE=/path/to/python
```

## Build and test

Use `gcc-debug` for normal development:

```bash
cmake --preset gcc-debug
cmake --build build/gcc-debug
ctest --test-dir build/gcc-debug --output-on-failure -E '^cpptest$'
```

For manual CLI use, install to a build-local prefix and run that executable; the build-tree CLI does not
have the bundled parameter sets:

```bash
cmake --install build/gcc-debug --prefix "$PWD/build/gcc-debug/_install"
build/gcc-debug/_install/bin/chargefw parameters
```

Gemmi registers its internal `cpptest` unconditionally even though that executable is excluded from
normal builds. Exclude it when invoking CTest directly.

After design is stable, build the affected target and run its focused regression test first. For example:

```bash
cmake --build build/gcc-debug --target test_methods
ctest --test-dir build/gcc-debug --output-on-failure -E '^cpptest$' -R '^test_methods$'
```

Native test executables use Snitch. To run a subset of one executable's test cases, pass a tag or a
test-name pattern; `--list-tags` and `--list-tests` show the available selections:

```bash
build/gcc-debug/tests/cpp/methods/test_methods '[sqeqp]'
```

Available presets are:

| Preset | Use |
| --- | --- |
| `gcc-debug` | Primary development and test loop |
| `clang-debug` | Cross-compiler validation |
| `gcc-release`, `clang-release` | Optimized and numerical validation |
| `clang-asan` | Memory and lifetime errors |
| `clang-ubsan` | Undefined behavior |
| `clang-tidy` | Static analysis during compilation |

All checked-in presets build the native library, CLI, tests, and Python bindings. Release presets also
enable IPO and host-native optimization; debug, sanitizer, and static-analysis presets disable
host-native optimization for reproducibility. Ignored `CMakeUserPresets.json` presets may inherit the
stable `gcc-debug` and `gcc-release` names to customize local build and installation directories.

Installation tests use the parent's build type, C and C++ compilers, compiler/linker flags, and IPO
settings for their nested builds. The custom-layout build also inherits ChargeFW's dependency, optimization,
sanitizer, warning, and developer-tool options; it changes the installation layout and builds only the
native library and CLI. The downstream consumer receives sanitizer requirements from the installed
library's exported target.

Configure another preset in the same way, then build and test its directory. `clang-tidy` runs during
its build and does not need a test run:

```bash
cmake --preset clang-tidy
cmake --build build/clang-tidy
```

Build sanitizer presets with low parallelism because they can use substantial memory:

```bash
cmake --preset clang-asan
cmake --build build/clang-asan --parallel 1
ctest --test-dir build/clang-asan --output-on-failure -E '^cpptest$'
```

On systemd-based Linux, an optional memory scope can protect an interactive workstation; adjust the
limits and parallelism for the machine:

```bash
systemd-run --user --scope \
    -p MemoryHigh=24G -p MemoryMax=32G -p MemorySwapMax=1G \
    cmake --build build/clang-asan --parallel 8
```

Choose validation by risk and qualification stage using [AGENTS.md](AGENTS.md#validation-cadence).
The full-suite and profile examples above document how to run selected checks, not a per-edit checklist.

## Python checks

CTest runs each `tests/python/test_<name>.py` file with pytest as `test_chargefw_python_<name>`, plus
the installed-package and mypy tests, with the correct build-tree `PYTHONPATH`. A focused Python test is
run like any other CTest test:

```bash
cmake --build build/gcc-debug --target chargefw_python
ctest --test-dir build/gcc-debug --output-on-failure -E '^cpptest$' \
    -R '^test_chargefw_python_calculation$'
```

To select individual pytest cases, run pytest directly against the build-tree package:

```bash
PYTHONPATH=build/gcc-debug/python python -m pytest tests/python/test_output.py -k failed_result
```

Ruff is configured in `pyproject.toml`. Run it on modified Python files:

```bash
python -m ruff check path/to/file.py
python -m ruff format --check path/to/file.py
```

## Formatting

Format modified C++ files with the repository configuration:

```bash
clang-format -i --style=file path/to/file.cpp
```

The pre-commit hook formats all tracked C and C++ files:

```bash
pre-commit run --all-files
```

## Continuous integration

`.github/workflows/ci.yml` runs on pull requests and pushes to `main`. It checks formatting with the
pre-commit hook using clang-format 22.1.8, runs Ruff over the repository, and builds and tests
`gcc-debug` with GCC 14 and Python 3.11 using the CTest command above. Update the pinned clang-format and
Ruff versions in the workflow together with the local tools.

The per-commit job uses the oldest supported toolchain on purpose: development usually happens with a
newer compiler and Python, which accept code that the oldest supported versions reject, so CI catches
that on every change. GCC 14 is the Debian 13 compiler. Pin the C and C++ compilers together: CMake
otherwise uses the runner's default C compiler, and link-time optimization fails when the two GCC
versions differ.

`.github/workflows/qualification.yml` runs weekly and on demand. It builds and tests the compatibility
containers (Debian 13 with the minimum supported dependency versions, Ubuntu 26.04, and Fedora 44 with
the newest packaged versions), builds and runs the CLI image, and tests the pinned dependency downloads
with `gcc-release`. It detects breakage that arrives without a commit, such as distribution package
updates; compiler, sanitizer, and clang-tidy checks stay in local release validation.

`.github/workflows/release.yml` is dispatched manually to build and publish a release; see
[RELEASING.md](RELEASING.md).

## Installation and containers

Use the [native source installation](docs/NATIVE.md#build-and-link), enabling `CHARGEFW_BUILD_CLI` when
the executable is required. It needs the dependency packages listed there, or
`CHARGEFW_FETCH_DEPENDENCIES=ON`, but not the Python and test dependencies used by the validation
presets.

The root `Dockerfile` builds the distributable CLI image:

```bash
docker build --tag chargefw:local .
```

For an image that will remain on the build machine, enable host-native optimization explicitly:

```bash
docker build \
    --build-arg CHARGEFW_ENABLE_NATIVE_OPTIMIZATIONS=ON \
    --tag chargefw:local-native .
```

The Dockerfiles under `docker/` build and test ChargeFW on Ubuntu 26.04, Debian 13, and Fedora 44. Run
these builds one at a time because they are resource intensive.

Python distribution builds and publication are documented in [RELEASING.md](RELEASING.md).
