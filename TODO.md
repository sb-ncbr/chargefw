# ChargeFW TODO

This file lists unfinished product, integration, scientific-qualification, distribution, and project
infrastructure work. Implemented behavior is documented under [`docs/`](docs/), and developer and release
procedures are in [DEVELOPMENT.md](DEVELOPMENT.md) and [RELEASING.md](RELEASING.md).

## 3.0.0 release

- [ ] Publish the Figshare item that archives each GitHub Release, and add its DOI to `CITATION.cff`
  and the README. Publishing the draft currently fails with HTTP 403. Check that the next release
  creates a new version of the same item. Complete the authors in `CITATION.cff` and add the paper
  reference once it exists.
- [ ] Reactivate the Read the Docs `stable` version at 3.0.0 and make it the default version.
- [ ] Check the documentation links on PyPI after the documentation reorganization. The default project
  page shows 0.1.7, whose description links to `docs/PYTHON.md`, `docs/FORMATS.md`, `docs/PROJECT.md`,
  and `docs/recipes/` on GitHub `main`; the 3.0.0a1 and 3.0.0a2 pages link to the old Read the Docs
  pages `PYTHON.html`, `FORMATS.html`, `PROJECT.html`, and `recipes/README.html`. All of them return
  404. Published descriptions cannot be edited. The 3.0.0 release replaces the default page; confirm
  that its links work and decide whether the remaining old pages need redirects or placeholder files.

## Distribution

- [ ] Add and qualify wheels for other platforms, architectures, and Python versions (Linux aarch64 and
  macOS arm64 first) as their native and optional Python dependencies permit, without significant extra
  machinery.
- [ ] Simplify Python source-build dependency staging and its loader handling once native external
  dependencies are qualified. The Python source build still downloads Gemmi and oneTBB and installs them
  next to the extension with the loader handling described in
  [DEVELOPMENT.md](DEVELOPMENT.md#dependencies). Decide, for example, between a static Gemmi and
  dependencies staged for wheel repair, and remove the custom loader rules once the replacement passes
  the source-build and wheel tests.

## Integrations

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

- [ ] Decide how to handle 2D coordinates. A 2D molfile or SD record, read by ChargeFW or converted from
  an RDKit conformer with `Is3D() == False`, becomes 3D coordinates with z = 0, and automatic selection
  then uses a geometry-dependent method on the flat drawing without a warning. Options: reject 2D
  coordinates for methods that use geometry, warn, or import the molecule without coordinates.

## Project infrastructure

- [ ] After 3.0.0, generate a Python API reference with `sphinx.ext.autodoc` from the published wheel
  installed on Read the Docs, without compiling ChargeFW there; the Python page remains the guide.
- [ ] After 3.0.0, add README badges for CI, documentation, PyPI version, and DOI.
- [ ] After 3.0.0, add `CONTRIBUTING.md`, linking to DEVELOPMENT.md and asking contributors to run the
  `gcc-debug` tests and pre-commit locally.
- [ ] After 3.0.0, enable Dependabot for GitHub Actions versions only.
