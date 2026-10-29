# ChargeFW TODO

This file lists unfinished product, integration, scientific-qualification, distribution, and project
infrastructure work. Implemented behavior is documented under [`docs/`](docs/), and developer and release
procedures are in [DEVELOPMENT.md](DEVELOPMENT.md) and [RELEASING.md](RELEASING.md).

## 3.0.0 release

- [ ] Reactivate the Read the Docs `stable` version at 3.0.0 and make it the default version.
- [ ] Drop the pre-release instructions at 3.0.0: remove `--pre` and the notes about it from both
  READMEs and the documentation, show the Docker image as plain `sbncbr/chargefw`, and remove the image
  tag update from `utils/release/01-prepare.sh` and RELEASING.md.
- [ ] Check the documentation links on PyPI after the documentation reorganization. The default project
  page shows 0.1.7, whose description links to `docs/PYTHON.md`, `docs/FORMATS.md`, `docs/PROJECT.md`,
  and `docs/recipes/` on GitHub `main`; the 3.0.0a1 and 3.0.0a2 pages link to the old Read the Docs
  pages `PYTHON.html`, `FORMATS.html`, `PROJECT.html`, and `recipes/README.html`. All of them return
  404. Published descriptions cannot be edited. The 3.0.0 release replaces the default page; confirm
  that its links work and decide whether the remaining old pages need redirects or placeholder files.

## Distribution

- [ ] Publish the Figshare item before the ChargeFW publication.

## Integrations

- [ ] ACC III migration: ChargeFW provides the reading options, method and parameter catalogs,
  applicability, calculation, and MOL2, PQR, and mmCIF output that ACC uses. Add an append mode to
  `chargefw.io.gemmi.attach_charges()` only if ACC needs to add charge sets to an existing file.

## Format correctness

- [ ] Decide explicit and hybrid bonds for `chargefw.io.gemmi.from_structure()`: Gemmi does not keep
  `_chem_comp_bond` or connection bond orders when writing a `Structure` as mmCIF, so ligand bond orders
  are lost. Options: refuse these strategies and default to `templates`, or document the loss.
- [ ] Decide whether to check generated mmCIF against the PDBx and charges dictionaries regularly, as a
  release step or a test. A manual check with `gemmi validate -d mmcif_pdbx_v50.dic -d
  mmcif_charges_v11.dic` passes; the optional parent checks (`-p`) report categories that minimal
  output omits, such as `_entity` and `_atom_type`.
- [ ] Leave atom names empty in the MOL/SDF reader, so generated output names are element plus atom
  number (`O1`, `O2`, `C3`) rather than the bare symbol followed by numbered duplicates (`O`, `O2`).
- [ ] After 3.0.0, compare the optional mmCIF model token exactly during strict attachment,
  distinguishing an absent column, `.`, `?`, and an explicit value.

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
