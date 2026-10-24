# ChargeFW CI and Documentation Plan (draft for review)

This is a proposal, not implemented behavior. After review, the accepted parts move into workflows,
[DEVELOPMENT.md](DEVELOPMENT.md), [RELEASING.md](RELEASING.md), and [TODO.md](TODO.md), and this file is
removed.

## Goals

- Protect `main` and pull requests with a fast, reliable check from a clean checkout.
- Keep expensive qualification (sanitizers, clang-tidy, release profiles, containers) off the per-commit
  path, following the [validation cadence](AGENTS.md#validation-cadence).
- Publish versioned user documentation built from the existing `docs/` Markdown.
- Replace token-based manual uploads with a workflow and trusted PyPI publishing, preserving the order
  and safeguards of [RELEASING.md](RELEASING.md): validate, build once, verify, TestPyPI, tag, approve,
  PyPI.
- Prepare the repository for public release and citation.

The repository is public, so GitHub-hosted runners are free and unmetered. The constraint is feedback time
on 2-vCPU runners, not cost.

## Overview

| Workflow | Trigger | Content | Expected role |
| --- | --- | --- | --- |
| `ci.yml` | PR, push to `main` | Lint, `gcc-debug` build and full CTest | Required merge check |
| `docs` (Read the Docs) | PR, push, tag | Sphinx build with warnings as errors | PR preview, published docs |
| `qualification.yml` | Weekly, manual | Distribution containers, CLI image, pinned dependencies | Drift detection |
| `release.yml` | Manual, on a `main` commit | sdist, wheels, TestPyPI, tag, approval, PyPI | Publication |

Local validation remains the primary qualification for design-sensitive PRs; CI is the clean-environment
safety net.

## 1. Per-commit CI (`ci.yml`)

Status: implemented in `.github/workflows/ci.yml` (lint with pinned clang-format 22.1.8 and Ruff 0.17.2;
`gcc-debug` on Ubuntu 26.04 with `g++-14`, Python 3.11, and ccache). The first run passed (lint 12 s, `gcc-debug` 6 min without cache), confirming GCC 14. Branch protection follows
the first green run.

Triggers: `pull_request`, `push` to `main`, `workflow_dispatch`.

Jobs:

1. **lint** (seconds): `pre-commit run --all-files` (clang-format) with a pinned `clang-format` version,
   and `ruff check` / `ruff format --check` on Python sources.
2. **gcc-debug** (main cost):
   - Ubuntu 26.04 runner with GCC 14 (the oldest GCC currently qualified by the Debian 13 container).
   - `actions/setup-python` with Python 3.11 (the minimum supported version), then the Python test
     dependencies from [DEVELOPMENT.md](DEVELOPMENT.md#requirements): NumPy, Nanobind, Gemmi 0.7.4,
     pytest, jsonschema, mypy. RDKit is not installed; its tests skip.
   - `cmake --preset gcc-debug -DCHARGEFW_PYTHON_EXECUTABLE=$(which python)`, build, and
     `ctest -E '^cpptest$'`. This covers native, CLI, Python, mypy, installed-package, relocation, and
     downstream-consumer tests.
   - ccache keyed on compiler and preset.
   - Dependencies are Ubuntu packages (including Gemmi 0.7.4) with `CHARGEFW_FETCH_DEPENDENCIES=OFF`;
     only the test framework is downloaded (see [section 6](#6-system-dependencies)).

Cost controls:

- `concurrency` with `cancel-in-progress` per branch.
- No workflow-level path filters: a required check whose workflow is skipped stays pending and blocks the
  PR. `ci.yml` always runs; documentation-only changes pay one cached build, which is acceptable at the
  expected PR volume. If that becomes slow, add a `changes` job and skip build steps inside the required
  jobs, which then still report success.
- Recipes are not run in CI; the release smoke check runs one.
- Minimal `permissions: contents: read`.

Branch protection on `main`: require `lint` and `gcc-debug`; allow maintainers to merge after local
qualification when CI is unavailable.

## 2. Documentation (Sphinx on Read the Docs)

Status: the basic site is published at <https://chargefw.readthedocs.io/> (`latest` from `main`,
default version). Implemented: Sphinx with MyST-Parser and Furo, `docs/index.md` as the landing page,
`docs/conf.py`, pinned `docs/requirements.txt`, `.readthedocs.yaml` with `fail_on_warning: true`, and
links from `pyproject.toml`, the PyPI README, and the GitHub README.

Remaining Read the Docs settings:

- `stable` is deactivated and the `v0.1.x` tags stay inactive; they predate the documentation
  configuration. Reactivate `stable` at 3.0.0 and make it the default version then.
- Add an automation rule "Activate version" for tags matching `^v\d+\.\d+\.\d+` (final and pre-release
  tags from 3.0.0a1 on). Versions are tied to Git tags, not GitHub Releases; the release workflow
  creates the tag after TestPyPI verification, so a tagged version's docs exist before its PyPI
  approval. An abandoned release keeps its docs version, which is deactivated manually.
- Enable pull-request preview builds.
- Add a second maintainer.

Plan:

- Sphinx with MyST-Parser so the existing `docs/*.md` remain the sources and continue to render on
  GitHub.
- Recipes are included verbatim with `literalinclude` on a small page per recipe or one recipes page;
  they are currently linked as downloads.
- Links to repository files outside `docs/` (README, TODO, sources) are absolute GitHub URLs so the
  Sphinx build is warning-free.
- Read the Docs builds `latest` from `main`, `stable` and per-version docs from tags, and preview builds
  for PRs. No GitHub workflow is needed for deployment. `stable` tracks the greatest non-pre-release tag,
  so `v3.0.0a1` appears only as its own version. New tags are created inactive; an automation rule
  activates tags matching `v*`.
- The Read the Docs Community plan limits builds to 15 minutes and 7 GB; the Markdown build and the
  phase 2 wheel installation fit well within it, which is why the extension is never compiled there.
- Phase 2 (after 3.0.0): Python API reference generated with `sphinx.ext.autodoc` by installing the
  published wheel on Read the Docs, instead of compiling the extension there. The hand-written
  [docs/PYTHON.md](docs/PYTHON.md) remains the guide; generated pages are the reference.
- The C++ API stays hand-documented in [docs/NATIVE.md](docs/NATIVE.md); no Doxygen.

The documentation deliverable for 3.0.0 is complete when:

- Read the Docs builds `latest` from `main`, PR previews, and `stable` from the 3.0.0 tag, all with
  warnings as errors;
- every page under `docs/`, including all recipes, is reachable from the site navigation;
- `pyproject.toml` `Documentation` and the README link to the site.

## 3. Qualification (`qualification.yml`)

Status: implemented in `.github/workflows/qualification.yml`; first green run for `ca3c314`.

Purpose: detect breakage that arrives without a commit (distribution package updates, base images,
unpinned Python packages, unavailable pinned downloads) and keep evidence for the documented dependency
versions. Compiler, sanitizer, clang-tidy, and release-profile checks are not repeated here; the
maintainer runs them locally before a release (see Decisions).

Triggers: weekly `schedule` and `workflow_dispatch`. It does not block merges; failures are reported
through GitHub notifications.

| Job | What it builds | Notes |
| --- | --- | --- |
| containers | `docker/*` | Debian 13 (minimum dependency versions), Ubuntu 26.04, and Fedora 44 (newest packaged versions), each with the full CTest |
| cli-image | `Dockerfile` | Builds the distributable CLI image and runs it |
| pinned-dependencies | `gcc-release` | Downloaded pins, which are newer than any distribution set (see [section 6](#6-system-dependencies)) |

Keep this workflow small; add a job only for a risk that neither `ci.yml` nor local validation covers.

## 4. Release (`release.yml`)

Versions: the numeric `project(VERSION)` plus `CHARGEFW_VERSION_SUFFIX` form the full PEP 440 version
used for metadata, runtime, output files, scripts, and tags (implemented; see
[RELEASING.md](RELEASING.md)).

Trigger: `workflow_dispatch` on `main` with the release version as input. The version must match the
repository's full distribution version and must not exist as a tag or on either package index. Version
preparation stays a reviewed commit made with `01-prepare.sh`.

Permissions are declared per job; the workflow default is `contents: read`. `testpypi` and `pypi` have
`id-token: write` for trusted publishing; `tag` and `github-release` have `contents: write`.

Jobs, in the order of [RELEASING.md](RELEASING.md):

1. **qualify**: require a successful `ci.yml` run on the dispatched commit. The Linux profiles are
   qualified locally with `02-validate.sh` before dispatch. Wheels for optional platforms (see
   Decisions) are tested on their own runners in the build job.
2. **build**: one sdist with `python -m build --sdist`, then wheels from that sdist with
   `pypa/cibuildwheel` using the existing `[tool.cibuildwheel]` configuration (CPython 3.11 to 3.14,
   manylinux x86-64, including the wheel test command). Optional platforms are added when they work
   without significant extra machinery (see Decisions).
3. **verify**: `twine check --strict`; expected artifact set (one sdist, one wheel per configured
   platform and Python version); third-party license notices present in every wheel and, as its
   contents require, the sdist; SHA-256 manifest; clean-environment install and
   smoke calculation (logic from `04-verify-local.sh`). The artifacts and manifest are uploaded once as
   a workflow artifact with long retention; later jobs download them and check the manifest.
4. **testpypi**: `pypa/gh-action-pypi-publish` with trusted publishing in a `testpypi` environment, with
   `skip-existing` only after hash comparison; then compare remote hashes, install from TestPyPI, and
   run the smoke check.
5. **tag**: check that production PyPI has no conflicting filenames, then create and push the annotated
   tag `vVERSION` on the dispatched commit.
6. **pypi**: `pypi` environment with required reviewer approval; publishes the same files, compares
   remote hashes, and verifies the PyPI install.
7. **github-release**: creates the GitHub Release from the tag with the artifacts, checksum manifest, and
   CHANGELOG section attached (pre-releases marked as such). Figshare archives it.
8. **docker-image**: builds the root `Dockerfile` from the tagged commit and pushes it to Docker Hub as
   `sb-ncbr/chargefw`, in a `dockerhub` environment holding the access token. Tags: `VERSION` for every
   release; `MAJOR.MINOR` and `latest` only for final releases. Runs after **pypi**, does not block the
   release, and is rerun on its own if the push fails.

Recovery: each job is a resumable phase. A job failing on an environmental, network, or index-propagation
problem is rerun with "Re-run failed jobs", which reuses the stored artifacts; build and verify are never
rerun after an upload. Upload jobs compare remote hashes before uploading, so a partial upload is safe to
resume. If source, version, or artifacts must change after any upload, the run is abandoned and the
release restarts with a new version. Tags are never moved or deleted.

Transition:

1. Release `3.0.0a1` with the scripts, which accept pre-release versions, to qualify the version
   model, Read the Docs tag activation, and Figshare on the established upload path.
2. Configure trusted publishers for `release.yml` on PyPI and TestPyPI, the `testpypi` and `pypi`
   GitHub environments, and the `sb-ncbr` Docker Hub organization with a `dockerhub` environment
   holding its access token.
3. Release `3.0.0a2` with the workflow; the scripts are the fallback until that release completes.
4. In one change after the successful release: rewrite [RELEASING.md](RELEASING.md) for the workflow,
   keep `01-prepare.sh` and `02-validate.sh` (local qualification), remove the build and publication
   scripts and stored upload tokens, and remove the corresponding [TODO.md](TODO.md) item.

The CLI container image and native installation license-notice checks are part of the license TODO item
and are verified by `qualification.yml` (`docker build .`) and the existing installation tests.

## 5. Public-release and citation items

Required for 3.0.0:

- `CITATION.cff` with authors, ORCIDs, and the paper reference once available.
- Figshare GitHub integration with auto-sync, so every GitHub Release becomes a new version of one
  public Figshare item with a versioned DOI (the base DOI cited in the paper). Status: a private draft
  item was created with the GitHub connector and auto-sync. The repository has no GitHub Releases yet,
  so the draft holds a `main` snapshot; it is published only once it contains the first 3.0.0
  pre-release, and the second release must be checked to create a new item version.

After 3.0.0:

- README badges: CI, documentation, PyPI version, DOI.
- `CONTRIBUTING.md` linking to DEVELOPMENT.md and stating that contributors run `gcc-debug` tests and
  pre-commit locally.
- Dependabot for GitHub Actions versions only.

## 6. System dependencies

> **Status: phase 1 implemented and accepted (see [Incremental steps](#incremental-steps)); the first
> `ci.yml` run with `CHARGEFW_FETCH_DEPENDENCIES=OFF` passed on GitHub for `11395c0`.** Phase 2 is not
> designed and needs its own review.

Goal: move closer to conventional, portable packaging without adding significant complexity. Linux
remains the required platform (see Decisions); this section must not introduce platform-specific
machinery, and it should remove some.

### Problem

[cmake/dependencies.cmake](cmake/dependencies.cmake) fetches pinned CLI11, nlohmann/json, Eigen,
nanoflann, oneTBB, and Gemmi unless `CHARGEFW_USE_SYSTEM_DEPENDENCIES=ON` (default `OFF`) finds a
compatible package; that option means "try system packages, then fetch". Fetched Gemmi and oneTBB are
shared libraries installed with ChargeFW, which requires ELF-specific loader handling so that the
private copies win over same-named system libraries (Gemmi's SONAME is unversioned):

- `$ORIGIN` RPATHs with `--disable-new-dtags` (DT_RPATH, ahead of `LD_LIBRARY_PATH`) in
  [cmake/install_rules.cmake](cmake/install_rules.cmake) and
  [python/CMakeLists.txt](python/CMakeLists.txt);
- native install rules for fetched dependency packages (`CHARGEFW_TBB_FETCHED`,
  `CHARGEFW_GEMMI_FETCHED`) and the bundled public json and Gemmi development packages;
- Python install rules that place the library and the fetched `gemmi_cpp` and `tbb` beside the
  extension ([python/CMakeLists.txt](python/CMakeLists.txt), the `gemmi_cpp` and `tbb` install
  blocks). This is what makes an ordinary `pip install` from the sdist work, because that path never
  runs wheel repair;
- tests that qualify this layout, including fetched source directories reused by
  [test_chargefw_custom_install_layout.cmake](tests/cpp/test_chargefw_custom_install_layout.cmake).

The installed package declares its dependency requirements separately from the build:
[chargefwConfig.cmake.in](cmake/chargefwConfig.cmake.in) hard-codes `nlohmann_json 3.12` and
`gemmi 0.7.4`.

### Evidence

Packaged versions (2026-10-24) against the current minimums:

| | CLI11 ≥2.7.2 | json ≥3.12 | Eigen ≥5.0 | nanoflann ≥1.14 | TBB ≥2023.1 | Gemmi ≥0.7.4 |
| --- | --- | --- | --- | --- | --- | --- |
| Ubuntu 26.04 | 2.6.1 | 3.11.3¹ | 3.4.0 | 1.9.0 | 2022.3 | 0.7.4 |
| Debian 13 | 2.4.1 | 3.11.3 | 3.4.0 | 1.7.1 | 2022.1 | — |
| Fedora 44 | 2.7.2 | 3.12.0 | 5.0.1 | 1.9.0 | 2022.3 | — |
| Current pins | 2.7.2 | 3.12.0 | 5.0.1 | 2.0.1 | 2023.1.0 | 0.7.4 |

¹ Ubuntu's package is versioned 3.12.0 but contains 3.11.3.

No distribution package for snitch exists in any of the three.

Experiment: Debian 13 container, `CHARGEFW_USE_SYSTEM_DEPENDENCIES=ON` with minimums lowered to the
Debian versions, Gemmi and snitch fetched, Python disabled, Debug build.

- All system packages were found.
- One compile failure: [json_output.cpp](src/adapters/native/json_output.cpp) assigns
  `std::optional<double>` (`execution_radius`, `ExecutionPolicy::radius()`) directly to JSON, which
  requires json 3.12. A local helper converting to a number or `null` fixes it (step 0 below).
- With that fix, 31 of 32 CTest tests pass, including the native method tests with Eigen 3.4 and
  nanoflann 1.7.
- `test_chargefw_cmake_consumer` fails because it requires the nlohmann/json package inside the
  installation, one of several bundled-layout assumptions.

The Debian 13 set is therefore a candidate minimum, not a qualified contract. Not yet tested: the
Python package, Clang, the release and sanitizer profiles, and Ubuntu 26.04 and Fedora 44.

### Contracts

Three build paths, each with its own contract:

1. **Native build and installation** (CLI, C++ library, containers, CI): every dependency is an
   external package found with `find_package`. Nothing is downloaded unless the user explicitly permits
   it. ChargeFW installs only its own files; whoever provides the dependencies also provides their
   runtime libraries.
2. **Python source build** (`pip install` from the sdist or a checkout): unchanged in phase 1. It
   downloads the pinned dependencies and installs `gemmi_cpp` and `tbb` beside the extension with the
   existing loader handling.
3. **Release wheels**: unchanged in phase 1. Built from the sdist by `cibuildwheel` and repaired as in
   [RELEASING.md](RELEASING.md).

Installed SDK compatibility is a separate contract from source-build minimums. ChargeFW's public
headers expose nlohmann/json and Gemmi types and link them publicly, so an installed ChargeFW requires
the same json and Gemmi versions it was built against. TBB, Eigen, nanoflann, and CLI11 are private:
consumers need the TBB runtime library, not its CMake package.

### Phase 1: native builds (implemented)

1. **Download permission.** Replace `CHARGEFW_USE_SYSTEM_DEPENDENCIES` with
   `CHARGEFW_FETCH_DEPENDENCIES`: `OFF` by default for native
   builds, so a missing or too-old dependency is a configuration error; `ON` in the Python build
   configuration in [pyproject.toml](pyproject.toml), preserving contracts 2 and 3. With it `ON`, only
   the pinned, checksum-verified downloads are used, as Python builds do today; there is no mixed
   "system, then download" mode. snitch remains a test-only download.
2. **Minimums.** Lower the build minimums to the Debian 13 set (CLI11 2.4, nlohmann/json 3.11.3,
   Eigen 3.4, nanoflann 1.7, TBB 2022.1; Gemmi stays 0.7.4), with the json 3.11 `std::optional` fix. The
   set becomes the documented contract once the acceptance criteria pass.
3. **Installed package.** Generate `find_dependency(... EXACT)` requirements for json and Gemmi from
   the versions the build actually used, so build and installed requirements cannot diverge. Do not add
   `find_dependency(TBB)` while TBB is linked privately into a shared library.
4. **Native installation.** With downloads disabled, install no dependency files and no dependency
   development packages. The fetched-dependency install rules and loader handling remain only for
   builds with downloads enabled (contract 2) until phase 2.
5. **Containers and CI.** The root `Dockerfile`, `docker/*`, and `ci.yml` install distribution packages
   and build with downloads disabled. Debian 13 and Fedora 44 build the pinned Gemmi into `/usr/local`
   in an earlier image stage. The CLI image's runtime stage installs the distribution's runtime
   libraries instead of copying bundled ones.

Acceptance criteria:

- full `gcc-debug` CTest with downloads disabled on Ubuntu 26.04 (CI) and in the Debian 13 and
  Fedora 44 containers;
- the minimum, distribution, and current-pin sets each build and pass the native tests; the
  current-pin set is one job with downloads forced, complementing the distribution jobs rather than
  testing every combination;
- a downstream consumer against an installation whose dependencies are in a non-default prefix,
  checking that discovery reaches the same json and Gemmi packages the build used;
- the custom-layout installation test passes with `OFF` (it reuses fetched source directories only
  when downloading, so it was left unchanged);
- the Python source build and release wheels pass unchanged, including a clean sdist installation
  verified after its build directory is removed.

### Incremental steps

Phase 1 lands as small commits, each reviewable and validated on its own. Behavior for users changes
only at step 6. Step 4 precedes step 5 because the containers and CI run the consumer test.

| Step | Change | Status |
| --- | --- | --- |
| 0 | Write optional execution radii to JSON explicitly (json 3.11 compatibility) | Committed (`334b343`) |
| 1 | `URL_HASH SHA256` for every pinned download, including snitch | Committed (`81718cd`) |
| 2 | Lower `find_package` minimums to the Debian 13 set (Eigen: 5, then 3.4, because each Eigen version file accepts one major version); checked in Debian 13 and Fedora 44 containers | Committed (`f425833`) |
| 3 | Add `CHARGEFW_FETCH_DEPENDENCIES` (default `ON`: pins only; `OFF`: installed packages only), replacing `CHARGEFW_USE_SYSTEM_DEPENDENCIES`; the Ubuntu compatibility image temporarily uses the pins | Committed (`b5e7e20`) |
| 4a | Consumer test uses the json and Gemmi packages the build found when they are not installed with ChargeFW (bundled json still required with `ON`); Ubuntu 26.04 with `OFF` passes 32/32 | Committed (`423a454`) |
| 4b | Installed package requires the exact json and Gemmi versions the build used (pin versions when downloading); Ubuntu 26.04 with `OFF` generates json 3.11.3 and Gemmi 0.7.4 and passes 32/32 | Committed (`2207e66`) |
| 4c | `INSTALL_RPATH_USE_LINK_PATH` on the library and CLI, so installed binaries find dependency libraries outside the system library paths (for example Gemmi in `/opt/deps`); project-built libraries stay `$ORIGIN`-relative. Ubuntu 26.04 with `OFF` passes 32/32 with Gemmi in `/opt/deps` and with packages in system paths | Committed (`df8223f`) |
| 4d | `INSTALL_RPATH_USE_LINK_PATH` on the Python extension as well: it links Gemmi directly, so a Python build against Gemmi outside the loader's search paths (Fedora's `/usr/local/lib64`) failed to import. Installed RPATH with downloads unchanged (`$ORIGIN`) | Committed (`651cdc4`) |
| 5a | Ubuntu 26.04: root `Dockerfile`, `docker/Dockerfile.ubuntu-26.04`, and `ci.yml` build with `OFF` and distribution packages. Ubuntu also needs `gemmi` and `zlib1g-dev`: `libgemmi-dev` ships a static `libgemmi_cpp.a` whose CMake targets reference the `gemmi` program. The CLI image's runtime stage adds `libtbb12`; the compatibility image passes 41/41 | Committed (`6c9f0c0`) |
| 5b | Debian 13 and Fedora 44 images build with `OFF`, with the pinned, checksum-verified Gemmi built into `/usr/local` before ChargeFW; both pass 41/41 (needs 4d) | Committed (`375db25`) |
| 6 | Default `OFF`; `ON` in the shared developer presets and in `pyproject.toml`; NATIVE.md lists the supported versions and installed-package requirements, DEVELOPMENT.md the preset behavior. Host `gcc-debug` 41/41 | Committed (`5c1c1f4`) |

Out of scope until step 6 is done: macOS and Windows experiments, phase 2, and removing loader handling.

### Phase 2: Python staging (later, separate change)

Only after phase 1 is qualified: decide how Python source builds stage Gemmi and TBB with less loader
handling (for example a static fetched Gemmi, or dependencies staged for wheel repair), and remove the
remaining custom loader rules if the replacement is qualified for contracts 2 and 3. Static linking of
ChargeFW into the extension is deferred; it does not remove the shared Gemmi and TBB dependencies and
adds another linkage configuration.

### Licenses

Discovery strategy does not change license obligations for binaries ChargeFW distributes. Header-only
Eigen, nanoflann, CLI11, and nlohmann/json contribute code to the CLI and library binaries, so native
binary distributions, the CLI image, and the wheels remain in the license audit. The sdist's notices
follow its inspected contents.

### Consequences for the plan

- `ci.yml`: Ubuntu 26.04 packages are the per-commit dependency set; the fetched-dependency cache is
  dropped for the native job.
- `qualification.yml`: containers cover the minimum (Debian 13) and distribution sets (Ubuntu 26.04,
  Fedora 44); one job covers the current pins with downloads forced.
- `release.yml`: unchanged by phase 1.
- Documentation: the dependency minimums and the installed-package requirements become user-facing
  contracts in [docs/NATIVE.md](docs/NATIVE.md), and the build requirements in
  [DEVELOPMENT.md](DEVELOPMENT.md).

### Decisions

- Download permission: `CHARGEFW_FETCH_DEPENDENCIES`, `OFF` for native builds and `ON` for Python builds.
- Minimum set: the Debian 13 versions, documented after qualification.
- Installed requirements: the exact json and Gemmi versions the build used.

## Rollout order

1. `ci.yml` and branch protection.
2. Sphinx documentation and Read the Docs project (basic site published; see section 2 for the
   remaining settings and recipe pages).
3. Phase 1 of [section 6](#6-system-dependencies): native dependency
   discovery, installed requirements, the json 3.11 fix, `ci.yml`, and the Dockerfiles; then the
   third-party license-notice item.
4. `qualification.yml` with its weekly schedule and the dependency-version containers (implemented).
5. `3.0.0a1` with the scripts (released); then trusted publishers, GitHub environments, and
   `release.yml`, with `3.0.0a2` via the workflow, following the transition in section 4.
6. Optional: Linux aarch64 and macOS arm64 wheels in `release.yml`, with platform-specific tests and
   artifact validation, if they fit without significant extra machinery. Not a 3.0.0 gate.
7. `CITATION.cff`, Figshare, the versioning policy, and the CLI image before 3.0.0; badges, contributor
   guidance, and autodoc after.

Section 6 phase 1 should be
implemented before `qualification.yml` and `release.yml` are written.

## Decisions

- The name stays "ChargeFW" everywhere; the version continues ChargeFW2's lineage. `chargefw` exists on
  PyPI (0.1, 0.1.2–0.1.7) and TestPyPI (0.1.1–0.1.7). Releases before 3.0.0 are PEP 440 pre-releases
  (`3.0.0a1`, `3.0.0b1`, `3.0.0rc1`); the compatibility policy after 3.0.0 is documented in
  [docs/PROJECT.md](docs/PROJECT.md) (TODO).
- Documentation is hosted on Read the Docs (Community plan, with its single unobtrusive ad), for
  versioned docs, PR previews, and no deployment workflow.
- Releases are dispatched manually; the workflow creates the tag after TestPyPI verification, as in the
  current procedure.
- Release tags are `vVERSION` (for example `v3.0.0a1`), as created by `08-tag.sh`.
- Sanitizers, clang-tidy, and the other Linux profiles are qualified locally by the maintainer before a
  release. CI release jobs cover only what cannot be run locally.
- Linux x86-64 is the required platform for 3.0.0 (wheels, native SDK, CLI). Other platforms (Linux
  aarch64, macOS, Windows) are welcome when they fit without significant extra complexity; they are not
  release gates. Once a platform's wheels are published, its CI build and tests gate the releases that
  publish them. Wheel platforms and native SDK/CLI platforms are qualified separately. Known portability
  points: library location through `dladdr` (`src/parameters/io/parameter_data_paths.cpp`) and peak
  memory through `getrusage` (`apps/chargefw/output.cpp`, which assumes Linux's kilobyte `ru_maxrss`;
  macOS reports bytes).
- The per-commit baseline is the oldest supported toolchain, GCC 14 and Python 3.11, with the C and C++
  compilers pinned together; newer compilers are covered locally and by the Fedora 44 container. The
  rationale is in [DEVELOPMENT.md](DEVELOPMENT.md#continuous-integration).
- The CLI image is published to Docker Hub only, as `sb-ncbr/chargefw`, the single image name in the
  documentation. It is built for linux/amd64; linux/arm64 follows the optional-platform rule above.
