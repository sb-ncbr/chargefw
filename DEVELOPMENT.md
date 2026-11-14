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
Test builds use an installed Snitch 1.3.2 when CMake finds one and otherwise download it, even with
`CHARGEFW_FETCH_DEPENDENCIES=OFF`.

The checked-in presets also build the Python bindings. Their Python interpreter needs Python 3.12 or
newer, development headers, NumPy 2.0 or newer, Nanobind 2.11 or newer, Gemmi 0.7.4, SciPy 1.18 or
newer for the fitting tests, and the test tools pytest, jsonschema, and mypy. RDKit is optional; its
tests are skipped when it is absent. The
presets use `/usr/bin/python3`; select another interpreter explicitly when needed:

```bash
cmake --preset gcc-debug \
    -DCHARGEFW_PYTHON_EXECUTABLE=/path/to/python
```

### Dependencies

ChargeFW has three build paths, each with its own dependency rules:

- **Native builds** (the library, command line, containers, and CI) use installed packages found with
  `find_package` and download nothing unless `CHARGEFW_FETCH_DEPENDENCIES=ON`. Without downloads, the
  installation contains only ChargeFW's own files; with them, it also contains the Gemmi and oneTBB
  runtime libraries, but not third-party development packages.
- **Python source builds** (`pip install` from the sdist or a checkout) always download the pinned,
  checksum-verified sources, as set in `pyproject.toml`. They link Gemmi statically into the extension
  and the packaged command line, because its library name carries no version and could meet another
  copy in the same process. oneTBB stays a shared library, installed once next to them under its
  interface version (`libtbb.so.12`); they find it through relative runtime search paths (`$ORIGIN` on
  Linux, `@loader_path` on macOS), and Windows loads DLLs from the extension's directory. On Linux,
  `--disable-new-dtags` makes this copy win over a same-named library from `LD_LIBRARY_PATH`. This path
  never runs wheel repair.
- **Release wheels** are built from the sdist by `cibuildwheel` and repaired as described in
  [RELEASING.md](RELEASING.md).

Third-party dependencies are private to ChargeFW. The JSON and mmCIF readers keep their parser state
behind private implementation pointers, so installed headers and CMake targets need no third-party
development packages. Internal targets that use a dependency directly link it explicitly. Consumers
still need any shared dependency libraries linked by ChargeFW at runtime. Native binaries, the CLI
image, and the wheels contain dependency code, so the license notices in [`licenses/`](licenses/)
ship with each of them.

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

All checked-in presets build the native library, CLI, tests, and Python bindings, and optimize for the
host CPU (`CHARGEFW_ENABLE_NATIVE_OPTIMIZATIONS=ON`), because they build and run on the same machine
and test speed matters. Release presets also enable IPO. Distributed builds use the x86-64-v3 default
instead, and the wheel workflow tests them. Ignored `CMakeUserPresets.json` presets may inherit the
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
`gcc-debug` with GCC 14, Python 3.12, and RDKit 2024.9 using the CTest command above. It also calls
`.github/workflows/wheels.yml`, which builds and tests the wheel for every published platform, so the
release workflow's requirement of a successful CI run covers all of them. Update the pinned clang-format
and Ruff versions in the workflow together with the local tools.

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

`.github/dependabot.yml` opens one pull request a month when the actions used by the workflows have new
major versions. Other dependencies are pinned and updated by hand.

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
these builds one at a time because they are resource intensive:

```bash
docker build --file docker/Dockerfile.ubuntu-26.04 .
docker build --file docker/Dockerfile.debian-13 .
docker build --file docker/Dockerfile.fedora-44 .
```

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
published: Linux x86-64 and aarch64, macOS arm64, and Windows x86-64 (see RELEASING.md). CI builds and
tests every published wheel, so each platform gates releases. Further platforms are welcome as wheels
when they need no significant extra machinery, and gate releases once their wheels are published. Intel
Macs are not planned: macOS 26 is their last release. Wheels and the native library are qualified
separately.

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

ChargeFW is library-first. File readers (`chargefw::io`) and the Python toolkit bridges translate
external representations into an owned, toolkit-neutral molecule model; method applicability,
parameter matching, execution selection, and calculation remain in the native library.

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

Imported records keep reader metadata beside the normalized molecule rather than inside
`core::Molecule`. The metadata owns only source correspondence, record-local import policy, diagnostics,
and a small connectivity summary; it does not retain source documents, parser objects, normalized-value
snapshots, or per-bond provenance.

The application-facing charge result owns those ordered imported records together with the toolkit-neutral
execution result. The writers consume that boundary directly, keeping source correspondence out of the
scientific calculation engine without requiring parallel caller-supplied structures or contexts.

### Parameterization

Fitting parameter values builds on the calculation without changing it. `chargefw::parameterization`
(`include/chargefw/parameterization/`, `src/parameterization/`) may use calculation, methods,
parameters, features, and core; none of them uses parameterization, and no public header outside
`chargefw/parameterization/` mentions fitting. Calculating a plan with the values of a compatible set,
`calculation::calculate(plan, candidate)`, stays in calculation: it is ordinary execution with other
values.

- **Reading a plan's assessment.** A plan shares the immutable assessment state
  (`src/calculation/assessment_state.h`) with its assessment and every copy.
  `src/calculation/assessment_access.h` gives the rest of the library read access, through a plan, to
  the assessed and prepared molecules, the fixed-ion partition, the plan's classifications, and the
  plan with candidate values. Keep it narrow, so that the public assessment header names no consumer.
- **Exact derivatives** of EEM, SQE, SQE+q0, and SQE+qp are private to the library
  (`src/parameterization/method_derivatives.h`). They are attached to the built-in methods from
  outside, matched by their concrete, final type; the methods know nothing of them. What they share
  with the methods is the system each method solves (`sqe_system.h`, `eem_system.h`) and each SQE
  variant's initial charges, written once for both. A change to a method's equations, such as its
  interaction function, needs the matching change to its derivatives; the central-difference test of
  every slot in `tests/cpp/parameterization/test_method_derivatives.cpp` catches one that is missed.
- **Derivatives per parameter-set value.** `normal_equations()` asks each molecule for the
  derivatives of all its charges with respect to the parameter-set values it uses, one column per
  value (`LinearizedCharges::jacobian()`, forward mode). That costs about one solve per value instead
  of one per atom, and each molecule's columns reduce at once to its part of Jᵀ W² J and Jᵀ W² r; the
  Jacobian of all charges is never formed.
- **Python.** The native fitting functions form the `_chargefw.parameterization` submodule
  (`python/src/parameterization.cpp`), apart from the calculation and parameters bindings. The package
  `chargefw.parameterization` keeps the optimizer, the profiles, and their decisions, such as which
  parameters of single-bonded atom types to hold; native code reports the facts those decisions rest
  on (`single_bond_entries()`). `_compression.py` turns the normal equations into the small
  least-squares model that SciPy's optimizer steps on.

The fitting profiles rest on these decisions:

- **Invariances are fixed, not found.** Some changes of values leave every charge unchanged: the same
  constant added to every EEM `A`, to every SQE electronegativity, or to every SQE+qp `q0`; EEM's `A`,
  `B`, and `kappa` scaled together; and the sign of an SQE width. Each profile fixes one slot per
  invariance, a default you can see and change, and
  `tests/cpp/parameterization/test_parameter_invariances.cpp` checks every declared invariance on the
  methods. Whatever else the data cannot determine shows in the diagnostics' rank; the fit reports it
  and does not repair it, so that its behavior stays predictable.
- **The SQE profiles guarantee a stable system.** Fitting the hardness as
  u = √π·hardness·width − 1 ≥ 0, with `width ≥ 0.05` and `kappa ≥ 0.01`, makes the interaction matrix
  the Coulomb energy of Gaussian charges plus a non-negative diagonal, so the split-charge matrix has
  its smallest eigenvalue at least the smallest `kappa`, for every molecule and geometry. The condition
  is sufficient, not necessary: published sets violate it and stay usable, but new fits meet it. EEM's
  point-charge interaction allows no such guarantee without a minimum distance between atoms.
- **An atom with one bond** exchanges charge only through that bond, so in the SQE methods its
  hardness and the `kappa` of its bond act as one sum. When those bonds' entries belong to its type
  alone and their `kappa` has no upper bound, the profile holds the hardness at u = 0 and lets the
  `kappa` take up the difference, which removes the direction without costing accuracy.
- **Validity limits stay in the methods; fitting bounds belong to the fit.** Methods reject only values
  where a formula breaks, such as EEM `kappa <= 0`. Bounds and other fit settings live in the profiles
  and `FitResult.bounds`, never in the fitted parameter file: how a set was fitted belongs in the
  publication that presents it.

Rejected along the way: derivatives from a reimplementation of each method in Python (it would drift
from the native one), replacing a plan's parameter-set pointer without checking compatibility, validity
ranges per parameter in method metadata (unknowable for most methods, and they would reject published
sets), choosing anchors for invariances automatically, L-BFGS-B (it moves a start outside the bounds
without notice), and a lower bound of 0 for strictly positive parameters (the optimizer's trial points
then reach values a method accepts but cannot use meaningfully).

### Departures from method publications

Where a method's publication does not cover an input, ChargeFW follows the method's construction and
says so in the method's "In ChargeFW" paragraph in `docs/guide/methods.md`. SQE+qp is the example:
Schindler et al. normalize the parameterized reference charges over all atoms of the molecule, and their
data sets had one connected component per molecule. Because SQE moves charge only along bonds, a single
normalization over several components would fix each component's total at a value set by the sizes and
parameters of the others. ChargeFW normalizes each component to its formal charge instead, as SQE+q0
does, which gives the published result for a single component.

### Issues and diagnostics

Each stage reports what it finds as its own issue type, with a kind enum and only the fields that
kind needs: `ImportIssue` from readers, `InputIssue` from assessment (and from execution for those that
concern only the executed method, such as `atoms_too_close`), `PrerequisiteIssue` and `ExecutionIssue`
from applicability, and `ChargeIssue` for calculated charges. `ResultDiagnostic` is the `io`
projection of all of them into one shape, with a string code and a severity, for result JSON and the
command line. Cautions never exclude a method; a finding that does is a `PrerequisiteIssue`.

A new kind follows these rules:

- **Bounded count.** One issue per kind and scope: record, molecule, or conformer, never one per atom.
- **Message for people.** It counts what was found and names one example, the most severe;
  `atom_index` is that atom.
- **Locations for programs.** An issue's `core::IssueLocation` lists every atom and bond it concerns in
  `atom_indices` and `bond_indices`, in index order, so that a viewer can select them. The message is
  never the only place a location appears. The lists hold at most the atoms and bonds of their scope,
  so they need no limit.
- **Stage.** Findings that do not depend on the method belong to assessment; those about the executed
  method to execution; those about charges exist only when charges do.
- **Python warning.** A kind also raises a Python warning (`InputWarning` or `ChargeWarning`, in
  `python/chargefw/_warnings.py`) when it can make the charges wrong and the user can act on it, by
  fixing the input or choosing another method. A kind that only notes data ChargeFW did not need,
  something the user requested, or a problem ChargeFW already refuses with a reason stays in its
  attribute.
