# ChargeFW TODO

This file lists unfinished product, integration, scientific-qualification, distribution, and project
infrastructure work. Implemented behavior is documented under [`docs/`](docs/), and developer and release
procedures are in [DEVELOPMENT.md](DEVELOPMENT.md) and [RELEASING.md](RELEASING.md).

## 3.0.0 release

The automation items under [Distribution](#distribution), the documentation site (with the completion
criteria in [CI.md](CI.md#2-documentation-sphinx-on-read-the-docs)), and `CITATION.cff` with Figshare
DOIs are also required for 3.0.0.

- [ ] Ship third-party license notices for bundled dependencies (shared Gemmi and oneTBB libraries;
  header-only Eigen, nanoflann, CLI11, nlohmann_json) in wheels, the sdist, native installations, and
  the CLI image.
- [ ] Build and qualify Linux aarch64 and macOS arm64 wheels.
- [ ] Document the versioning and compatibility policy in [docs/PROJECT.md](docs/PROJECT.md).
- [ ] Start `CHANGELOG.md` at the first 3.0.0 pre-release, stating the succession from ChargeFW2.

## Distribution

- [ ] Automate the manual `cibuildwheel` release process in [RELEASING.md](RELEASING.md), including
  artifact retention and publication through trusted PyPI environments.
- [ ] Automate repository validation, including formatting, GCC and Clang debug/release tests,
  AddressSanitizer and UndefinedBehaviorSanitizer tests, installed-package relocation, and the downstream
  CMake consumer.
- [ ] Add and qualify wheels for other platforms, architectures, and Python versions as their native and
  optional Python dependencies permit.

## Integrations

- [ ] Qualify the optional RDKit adapter against supported RDKit releases with real-toolkit conversion,
  mapping, property attachment, and SD serialization tests.
- [ ] Compare the Python API with the capabilities actually used by the ACC III backend and add only the
  missing explicit capabilities needed for migration. Document intentionally unsupported legacy
  behavior before replacing the existing backend.

## Format correctness

- [ ] Preserve the exact optional mmCIF model token during strict attachment, distinguishing an absent
  column, `.`, `?`, and an explicit value.
- [ ] Define and enforce how structural import handles Gemmi connections whose endpoints are in different
  symmetry images; the core graph cannot represent a symmetry-image endpoint.
- [ ] Make mmCIF charge annotation update stale or duplicate `_audit_conform` dictionary entries
  deterministically.
- [ ] Add optional dictionary-backed validation of generated mmCIF to focused or release validation.
- [ ] Define and test explicit-connectivity behavior for `chargefw.io.gemmi.from_structure()`, accounting
  for connectivity that Gemmi does not retain when converting a `Structure` through mmCIF.

## Scientific qualification

- [ ] Establish method-specific validation and documented tolerances before making broad accuracy claims
  for cutoff or cover execution. Current whole-molecule-radius and representative-structure checks verify
  implementation behavior, not a general reduced-mode accuracy envelope.
