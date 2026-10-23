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

GitHub-hosted runners are free and unmetered for public repositories. The constraint is feedback time on
2-vCPU runners, not cost. Before the repository is public, the private-repository quota (2,000
minutes/month) applies, so heavy workflows remain manual until then.

## Overview

| Workflow | Trigger | Content | Expected role |
| --- | --- | --- | --- |
| `ci.yml` | PR, push to `main` | Lint, `gcc-debug` build and full CTest | Required merge check |
| `docs` (Read the Docs) | PR, push, tag | Sphinx build with warnings as errors | PR preview, published docs |
| `qualification.yml` | Weekly, manual, release | Clang, release, sanitizers, clang-tidy | Drift detection, release gate |
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
   - ccache keyed on compiler and preset, and a cache of CMake-fetched dependencies.

Cost controls:

- `concurrency` with `cancel-in-progress` per branch.
- No workflow-level path filters: a required check whose workflow is skipped stays pending and blocks the
  PR. `ci.yml` always runs; documentation-only changes pay one cached build, which is acceptable at the
  expected PR volume. If that becomes slow, add a `changes` job and skip build steps inside the required
  jobs, which then still report success.
- Recipes are not run by `ci.yml`; recipe checks belong to `qualification.yml` (section 3).
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

Triggers: weekly `schedule`, `workflow_dispatch` (any branch, used for risky PRs), and `workflow_call`
from `release.yml`.

Matrix (independent jobs, `fail-fast: false`):

| Job | Preset | Notes |
| --- | --- | --- |
| clang-debug | `clang-debug` | Pinned LLVM from apt.llvm.org |
| gcc-release | `gcc-release` | Host-native optimization is acceptable on the test runner |
| rdkit | `gcc-debug` | Matrix over the oldest supported (`2024.9`) and latest RDKit; qualifies the [RDKit TODO](TODO.md#integrations) |
| clang-release | `clang-release` | |
| clang-asan | `clang-asan` | Build with `--parallel 2` |
| clang-ubsan | `clang-ubsan` | Build with `--parallel 2` |
| clang-tidy | `clang-tidy` | Build only |
| containers | `docker/*` | Matrix over Ubuntu 26.04, Debian 13, Fedora 44; also `docker build .` |
| recipes | installed package | Optional: runs each recipe on a test fixture, if recipes gain a fixture-based invocation |

The scheduled run reports failures through GitHub notifications; it does not block merges. These
checks remain runnable locally as documented in [DEVELOPMENT.md](DEVELOPMENT.md).

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

1. **qualify**: `ci.yml` jobs and `qualification.yml` on the dispatched commit (see open question 2).
2. **build**: one sdist with `python -m build --sdist`, then wheels from that sdist with
   `pypa/cibuildwheel` using the existing `[tool.cibuildwheel]` configuration (CPython 3.11 to 3.14,
   manylinux x86-64, including the wheel test command). Linux aarch64 and macOS arm64 runners are added
   as a separate milestone before `3.0.0b1` (see Rollout order).
3. **verify**: `twine check --strict`; expected artifact set (one sdist, one wheel per configured
   platform and Python version); third-party
   license notices present in the sdist and every wheel; SHA-256 manifest; clean-environment install and
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

Recovery: each job is a resumable phase. A job failing on an environmental, network, or index-propagation
problem is rerun with "Re-run failed jobs", which reuses the stored artifacts; build and verify are never
rerun after an upload. Upload jobs compare remote hashes before uploading, so a partial upload is safe to
resume. If source, version, or artifacts must change after any upload, the run is abandoned and the
release restarts with a new version. Tags are never moved or deleted.

Transition:

1. Release `3.0.0a1` with the scripts, which accept pre-release versions, to qualify the version
   model, Read the Docs tag activation, and Figshare on the established upload path.
2. Configure trusted publishers for `release.yml` on PyPI and TestPyPI, and the `testpypi` and `pypi`
   GitHub environments.
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

## Rollout order

1. `ci.yml` and branch protection.
2. Sphinx documentation and Read the Docs project (basic site published; see section 2 for the
   remaining settings and recipe pages).
3. `qualification.yml`, manual only while private; enable the schedule when public.
4. `3.0.0a1` with the scripts; then trusted publishers, GitHub environments, and `release.yml`, with
   `3.0.0a2` via the workflow, following the transition in section 4.
5. Linux aarch64 and macOS arm64 wheels in `release.yml`, with platform-specific tests and artifact
   validation, before `3.0.0b1`.
6. `CITATION.cff` and Figshare before 3.0.0; badges, contributor guidance, and autodoc after.

`ci.yml`, documentation, and `qualification.yml` do not depend on the open questions below; questions
1–3 must be settled before `release.yml` is written.

## Decisions

- The name stays "ChargeFW" everywhere; the version continues ChargeFW2's lineage. `chargefw` exists on
  PyPI (0.1, 0.1.2–0.1.7) and TestPyPI (0.1.1–0.1.7). Releases before 3.0.0 are PEP 440 pre-releases
  (`3.0.0a1`, `3.0.0b1`, `3.0.0rc1`); the compatibility policy after 3.0.0 is documented in
  [docs/PROJECT.md](docs/PROJECT.md) (TODO).
- Documentation is hosted on Read the Docs (Community plan, with its single unobtrusive ad), for
  versioned docs, PR previews, and no deployment workflow.
- Releases are dispatched manually; the workflow creates the tag after TestPyPI verification, as in the
  current procedure.

## Open questions

1. Tag format: keep `vVERSION` (current `08-tag.sh`)?
2. Should releases run `qualification.yml` in CI, or is the maintainer's local `02-validate.sh` run
   sufficient (shorter release workflow)?
3. Should the CLI container image be published (e.g. GHCR) by `release.yml`?
4. Is GCC 14 / Python 3.11 the right per-commit baseline, or should it track the newest toolchain?
