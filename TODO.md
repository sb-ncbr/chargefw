# ChargeFW TODO

This file lists unfinished product, integration, scientific-qualification, distribution, and project
infrastructure work. Implemented behavior is documented under [`docs/`](docs/), and developer and release
procedures are in [DEVELOPMENT.md](DEVELOPMENT.md) and [RELEASING.md](RELEASING.md).

## 3.0.0 release

- [ ] Reactivate the Read the Docs `stable` version at 3.0.0 and make it the default version. Do it
  after the release workflow has pushed the tag and before approving the PyPI upload, so that the
  `/en/stable/` links of the package description and the `$id` addresses of the JSON Schemas work
  from the start.
- [ ] Set the `Development Status` classifier in `pyproject.toml` to `5 - Production/Stable` at 3.0.0.
- [ ] Check the documentation links of the 3.0.0 page on PyPI, which replaces 0.1.7 as the default
  project page. The 3.0.0a1 and 3.0.0a2 pages link to the old Read the Docs pages `PYTHON.html`,
  `FORMATS.html`, `PROJECT.html`, and `recipes/README.html`, which return 404; decide whether they
  need redirects.

## Distribution

- [ ] When Gemmi releases a version that builds with nanobind 3 or has CPython 3.15 wheels
  ([project-gemmi/gemmi#471](https://github.com/project-gemmi/gemmi/issues/471)), move the
  Gemmi pins to it (`pyproject.toml`, `cmake/dependencies.cmake`, the CI and container installs, and
  the license notice) and remove the nanobind build constraint from `tests/python/build-constraints.txt`,
  `pyproject.toml`, `docker/Dockerfile.fedora-45`, the installation guide, and RELEASING.md.

- [ ] Publish the Figshare item before the ChargeFW publication, and add a badge with its DOI to
  `README.md`.

## Integrations

- [ ] ACC III migration: the development ACC runs on ChargeFW 3, and production still uses
  ChargeFW2. [ACC.md](ACC.md) has the working notes; remove it when production has moved.

## Parameterization

- [ ] Fit to several targets at once: electrostatic potentials at points, and the response of the
  charges to external charges, which needs a second plan with fixed sources. The targets share one
  parameter layout (`parameters::compatible()`); `fit()` takes a list of them and reports diagnostics
  per target, and `fit(evaluator, references)` stays the shorthand for one charge target. Planned for
  3.1.0: branch `multi-target` holds the difference normal equations, their Python binding, and
  explicit per-atom fixed charges in `fixed_ions`, which need a new attribute in `RequestedCalculation`
  and `CalculationSettings`.
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
