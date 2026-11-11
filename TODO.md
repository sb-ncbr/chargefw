# ChargeFW TODO

This file lists unfinished product, integration, scientific-qualification, distribution, and project
infrastructure work. Implemented behavior is documented under [`docs/`](docs/), and developer and release
procedures are in [DEVELOPMENT.md](DEVELOPMENT.md) and [RELEASING.md](RELEASING.md).

## 3.0.0 release

- [ ] Reactivate the Read the Docs `stable` version at 3.0.0 and make it the default version.
- [ ] Set the `Development Status` classifier in `pyproject.toml` to `5 - Production/Stable` at 3.0.0.
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

## Calculation

- [ ] Python: reuse the choices of a plan on other molecules. A plan would give its fully explicit
  settings (method, parameter set, options, parameter matching, execution mode and radius, fixed
  ions), which `calculate()` and `assess()` accept for new molecules; they apply exactly those
  choices and fail instead of falling back. Needs a design sketch and review first.

## Parameterization

- [ ] Fit to several targets at once: electrostatic potentials at points, and the response of the
  charges to external charges, which needs a second plan with fixed sources. The targets share one
  parameter layout (`parameters::compatible()`); `fit()` takes a list of them and reports diagnostics
  per target, and `fit(evaluator, references)` stays the shorthand for one charge target.
- [ ] Charged groups: find resonance groups and formal charges shared within them (carboxylate,
  guanidinium) from the bond graph into `features::TopologyFeatures`, independent of atom order and of
  the resonance form in the input; add an ionic-group atom classifier and a resonance bond
  classification; and give SQE with shared initial charges its own method ID, because sets fitted with
  sharing are valid only with it. Decide first which groups count.
- [ ] Profile `fit()` with many small molecules: a fit to 25,628 fragments of at most 257 atoms used
  about 6 of 16 cores. Measure where the time goes before changing anything.
- [ ] Investigate an atom classifier that sees beyond the first neighbors, such as the neighbors' bond
  orders and formal charges or the second shell, as a new classifier value. Hydroxyl oxygens of
  carboxylic acids and of alcohols share every first-neighbor type. Whether it is worth adding depends
  on the accuracy finer types can give, which is still to be measured.
- [ ] Investigate metal–ligand coordination as bonds, so that SQE-family methods can move charge
  between a bound metal and its ligands instead of holding the metal as a fixed ion. This is a
  chemistry policy: no coordination bonds are inferred automatically until the policy is decided.

## Project infrastructure

- [ ] After 3.0.0, generate a Python API reference with `sphinx.ext.autodoc` from the published wheel
  installed on Read the Docs, without compiling ChargeFW there; the Python page remains the guide.
- [ ] After 3.0.0, add README badges for CI, documentation, PyPI version, and DOI.
- [ ] After 3.0.0, add `CONTRIBUTING.md`, linking to DEVELOPMENT.md and asking contributors to run the
  `gcc-debug` tests and pre-commit locally.
- [ ] After 3.0.0, enable Dependabot for GitHub Actions versions only.
