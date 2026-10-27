# ChargeFW development

This is the practical guide to building and testing ChargeFW, with an [architecture](#architecture)
overview at the end. Repository policy and architecture rules are in [AGENTS.md](AGENTS.md); public
behavior is documented under [`docs/`](docs/).

## Requirements

Native development requires CMake 3.28 or newer, Ninja, and a C++23-capable GCC or Clang toolchain.
The checked-in presets set `CHARGEFW_FETCH_DEPENDENCIES=ON`, so CMake downloads pinned,
checksum-verified dependency sources. CI and the compatibility containers build against distribution
packages with `CHARGEFW_FETCH_DEPENDENCIES=OFF`; pass it to a preset to use installed packages locally.
The supported versions are listed in [the C++ page](docs/interfaces/cpp.md#build-and-link).

The checked-in presets also build the Python bindings. Their Python interpreter needs Python 3.12 or
newer, development headers, NumPy, Nanobind, Gemmi 0.7.4, and the test tools pytest, jsonschema, and
mypy. RDKit is optional; its tests are skipped when it is absent. The presets use `/usr/bin/python3`;
select another interpreter explicitly when needed:

```bash
cmake --preset gcc-debug \
    -DCHARGEFW_PYTHON_EXECUTABLE=/path/to/python
```

### Dependencies

ChargeFW has three build paths, each with its own dependency rules:

- **Native builds** (the library, command line, containers, and CI) use installed packages found with
  `find_package` and download nothing unless `CHARGEFW_FETCH_DEPENDENCIES=ON`. Without downloads, the
  installation contains only ChargeFW's own files; with them, it also contains the nlohmann/json and
  Gemmi development packages and the oneTBB runtime.
- **Python source builds** (`pip install` from the sdist or a checkout) always download the pinned,
  checksum-verified sources, as set in `pyproject.toml`, and install the fetched Gemmi and oneTBB
  libraries next to the extension, which finds them through relative runtime search paths (`$ORIGIN` on
  Linux, `@loader_path` on macOS); Windows loads DLLs from the extension's directory. On Linux,
  `--disable-new-dtags` makes these copies win over same-named libraries from `LD_LIBRARY_PATH`;
  Gemmi's library name carries no version. This path never runs wheel repair.
- **Release wheels** are built from the sdist by `cibuildwheel` and repaired as described in
  [RELEASING.md](RELEASING.md).

ChargeFW's public headers expose nlohmann/json and Gemmi types, so the installed CMake package requires
the exact json and Gemmi versions the build used; `cmake/chargefwConfig.cmake.in` generates those
requirements. oneTBB, Eigen, nanoflann, and CLI11 are private: consumers need only the oneTBB runtime
library. Native binaries, the CLI image, and the wheels all contain code of these dependencies, so their
license notices in [`licenses/`](licenses/) ship with each of them.

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

Build and install a local wheel from the source tree with:

```bash
uv build --quiet --wheel
uv pip install --link-mode=copy --reinstall dist/chargefw-*.whl
```

The wheel contains the private native extension, required shared libraries, bundled parameter JSON,
Python modules, type declarations, and the `py.typed` marker. Parameter discovery uses package resources
and does not depend on the current directory or an environment variable.

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
`gcc-debug` with GCC 14, Python 3.12, and RDKit 2024.9 using the CTest command above. Update the pinned
clang-format and Ruff versions in the workflow together with the local tools.

The per-commit job uses the oldest supported toolchain on purpose: development usually happens with a
newer compiler and Python, which accept code that the oldest supported versions reject, so CI catches
that on every change. GCC 14 is the Debian 13 compiler. Pin the C and C++ compilers together: CMake
otherwise uses the runner's default C compiler, and link-time optimization fails when the two GCC
versions differ.

`.github/workflows/qualification.yml` runs weekly and on demand. It builds and tests the compatibility
containers (Debian 13 with the minimum supported dependency versions, Ubuntu 26.04, and Fedora 44 with
the newest packaged versions), builds and runs the CLI image, and tests the pinned dependency downloads
with `gcc-release` and the newest RDKit. It detects breakage that arrives without a commit, such as
distribution package updates or a new RDKit release; compiler, sanitizer, and clang-tidy checks stay in
local release validation.

`.github/workflows/release.yml` is dispatched manually to build and publish a release; see
[RELEASING.md](RELEASING.md).

The `main` ruleset requires `lint` and `gcc-debug`; repository admins can bypass it, for example to
merge after local qualification when CI is unavailable. The workflows have no path filters on purpose:
a required check whose workflow is skipped stays pending and blocks the pull request. If
documentation-only changes ever make CI too slow, skip the build steps inside the required jobs instead,
so they still report success.

## Installation and containers

Use the [native source installation](docs/interfaces/cpp.md#build-and-link), enabling `CHARGEFW_BUILD_CLI` when
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

## Documentation

The user documentation in [`docs/`](docs/) is built with Sphinx, MyST-Parser, and the Furo theme and
published on Read the Docs at <https://chargefw.readthedocs.io/>. Follow the
[writing rules](AGENTS.md#writing-user-documentation) when editing it. Build it locally with warnings as
errors, as Read the Docs does; `-n` additionally checks every cross-reference:

```bash
python -m venv /tmp/docs-venv
/tmp/docs-venv/bin/pip install -r docs/requirements.txt
/tmp/docs-venv/bin/sphinx-build -W -n -b html docs /tmp/docs-html
```

`docs/_ext/parameter_sets.py` generates the parameter-set tables from `data/parameters` and the method
headers into the ignored `docs/_generated/` directory before each build. Read the Docs never compiles
ChargeFW, which keeps its builds within the 15-minute limit of the free plan.

Read the Docs builds `latest` from `main` and a preview for every pull request. An automation rule
activates every tag matching `^v\d+\.\d+\.\d+`, so each release, including pre-releases, gets its own
documentation version. `stable` follows the newest final release; it stays inactive until 3.0.0. A
version of an abandoned release is deactivated by hand.

## Platforms

The C++ library and its native installation, including the native command line, are supported on Linux.
The Python package, which includes its own `chargefw` command line, is supported wherever wheels are
published. Linux x86-64 is required for both. Other platforms, such as Linux aarch64, macOS, or Windows,
are welcome as wheels when they need no significant extra machinery, but they do not gate releases until
their wheels are published; from then on, their CI build and tests do. Wheels and the native library are
qualified separately.

Wheel builds set `CHARGEFW_BUILD_LIBRARY=OFF` in `pyproject.toml`: the extension and the packaged command
line link the library objects, so the shared library and its installation are not built. Native builds
on macOS work for development but are not supported. CMake refuses them on Windows, because the library
exports no symbols for a DLL.

Operating-system calls live in `src/platform/`. Each platform directory (`linux/`, `macos/`, `windows/`)
implements the headers there under the same file names, and `src/platform/CMakeLists.txt` is the only
place that selects one; the rest of the code calls `chargefw::platform` without platform conditionals.
Platform build settings are in `cmake/`: runtime search paths and the ELF-only `--disable-new-dtags` in
`runtime_paths.cmake`, and compiler-specific warning flags in `warnings.cmake`.

macOS builds use Xcode's compiler and Apple's libc++, whose runtime is part of macOS. Two C++ library
features are therefore avoided, each with a comment where it would be used: `std::ranges::iota`, which
the libc++ in Xcode 26.1 lacks, is written as `std::iota`, and floating-point `std::from_chars`, which
Apple's libc++ provides only from macOS 26, is replaced by `gemmi::fast_from_chars`. The code then needs
macOS 13.3, the first with floating-point `std::format`; wheel tags record only the major version, so
the macOS wheels target 14.0. Restore `std::ranges::iota` once the Xcode used for the wheels provides
it, and `std::from_chars` once the oldest supported macOS is 26. New C++ library features need the same
check against Apple's libc++.

## Architecture

ChargeFW is library-first. File readers and toolkit adapters translate external representations into an
owned, toolkit-neutral molecule model; method applicability, parameter matching, execution selection,
and calculation remain in the native library.

**Algorithms are reusable and concurrent.** Prepared assessments and concrete plans own the state they
need and can be executed repeatedly without global mutable calculation state.

```text
files, C++ values, or Python arrays
                 |
                 v
       toolkit-neutral molecules
       atoms + bonds + conformers
                 |
                 v
       prepared topology/geometry
                 |
        +--------+--------+
        |                 |
        v                 v
 bundled parameters   method registry
        |                 |
        +--------+--------+
                 v
       applicability assessment
                 |
                 v
       concrete execution plans
                 |
                 v
          full/cutoff/cover
                 |
                 v
      source-ordered charge sets
```

Application-facing calculation functions compose the owned assessment facade rather than exposing the
lower-level prepared-feature and classification interfaces. Callers need an explicit assessment only when
they want to inspect or reuse its concrete plans.

Imported records keep adapter metadata beside the normalized molecule rather than inside
`core::Molecule`. The metadata owns only source correspondence, record-local import policy, diagnostics,
and a small connectivity summary; it does not retain source documents, parser objects, normalized-value
snapshots, or per-bond provenance.

The application-facing charge result owns those ordered imported records together with the toolkit-neutral
execution result. Output adapters consume that boundary directly, keeping source correspondence out of the
scientific calculation engine without requiring parallel caller-supplied structures or contexts.
