# ChargeFW Repository Change Policy

## Required reading and document ownership

Read this file before changing the repository. Then read only the relevant sections of:

- [docs/PROJECT.md](docs/PROJECT.md): implemented architecture, philosophy, capabilities, and limits;
- [docs/FORMATS.md](docs/FORMATS.md): language-independent molecular input and charge output behavior;
- [docs/CLI.md](docs/CLI.md): implemented command-line behavior;
- [docs/NATIVE.md](docs/NATIVE.md): implemented public C++ API and installation;
- [docs/PYTHON.md](docs/PYTHON.md): implemented Python API and package behavior;
- [docs/PARAMETERS.md](docs/PARAMETERS.md): parameter-set JSON and classifier behavior;
- [docs/recipes/](docs/recipes/): executable end-user examples, not normative behavior;
- [DEVELOPMENT.md](DEVELOPMENT.md): executable local build, test, analysis, and container workflows;
- [RELEASING.md](RELEASING.md): Python distribution validation and publication procedure;
- [TODO.md](TODO.md): unfinished product, integration, qualification, and infrastructure work; and
- [README.md](README.md): concise project entry point.

Keep ownership clear. The files under `docs/` are user documentation and must describe implemented
behavior, not development plans or milestone history. Unfinished work belongs in `TODO.md`; remove it
when complete rather than retaining checked history. `README.md` should remain a concise user entry
point. `DEVELOPMENT.md` and `RELEASING.md` own executable contributor and maintainer procedures;
`AGENTS.md` owns repository policy and should link to those procedures rather than duplicate them.
Document cross-interface concepts, rationale, and limitations once in `docs/PROJECT.md`. Keep molecular
serialization behavior in `docs/FORMATS.md` and parameter-file behavior in `docs/PARAMETERS.md`.
Interface-specific documents should state only the syntax and behavior needed for that interface and link
to the owning document when additional explanation is useful; do not repeat the same explanatory text
across CLI, native, and Python documentation.

Describe user-facing behavior through its positive contract: state what users can rely on and how to use
it. Avoid documenting absent fields, unsupported cases, or missing capabilities as standalone negative
contracts; document a limitation only when it is necessary to use the supported behavior safely, and frame
it around the supported boundary or alternative.

## Repository map

```text
include/chargefw/   Public C++ API
src/                Native library implementation, organized by public domain
apps/chargefw/      Command-line application
python/chargefw/    Public Python package and private extension stubs
python/src/         Nanobind extension implementation
data/parameters/    Bundled parameter sets
tests/cpp/          Native tests and downstream CMake consumer
tests/cpp/support/  Shared native test fixtures and CLI test script helpers
tests/python/       Python package and adapter tests
tests/fixtures/     Molecular and parameter test data
cmake/              Dependencies, diagnostics, and installation rules
docs/               Implemented user documentation
licenses/           Third-party license notices shipped with binaries
DEVELOPMENT.md       Developer build and validation manual
RELEASING.md         Maintainer package release procedure
Dockerfile           User-facing CLI container image
docker/              Developer compatibility containers
```

`include/chargefw/adapters/{native,gemmi}/all.h` are application convenience headers. Library code
should include only the individual adapter headers it uses.

## Architectural rules

- Keep `core` small and toolkit-neutral. File formats, toolkit objects, preparation policy, and cached
  derived data do not belong in `core::Molecule`.
- Put cached topology/geometry data and spatial fragments in `features`.
- Parameter matching is immutable external data:
  `ParameterSet + Molecule + TopologyFeatures -> ParameterClassification -> ParameterView`.
- Methods are stateless algorithms. Request-specific options, geometry, targets, and parameters enter
  through `methods::CalculationInput`.
- Keep scientific applicability, execution availability, deterministic selection, and calculation
  distinct, while allowing the application facade to compose them.
- Keep execution policy separate from method options. Never silently change method, parameter set,
  classification, execution mode/radius, topology, protonation, or geometry.
- Preserve source atom order, molecule/conformer identity, and mappings at every result boundary.
- Missing scientific prerequisites are hard failures. Resource thresholds only guide automatic
  execution; explicit full execution may override them with a reported warning.
- Full calculations must remain faithful to their cited methods. Approximate execution must be
  capability-checked, explicit in provenance, and validated per method.
- Bindings, adapters, and applications translate representations into the native facade; they must not
  duplicate applicability, selection, or scientific policy.
- Avoid global mutable state. The library must remain suitable for concurrent and embedded use.

## Change discipline

Use existing context and inspect the affected code and contracts. Delegate substantial work with a concise,
phase-specific brief and one owner; reuse that worker for the same slice. Handle small local edits directly
and avoid duplicate exploration. Handoffs should report the diff, relevant evidence, and unresolved issues.

For design-sensitive API, ownership, scientific, or representation changes, first produce a small
representative code sketch and request user architectural review before full tests or propagation to
bindings and documentation. Label unqualified designs as drafts. A small experiment or compile check is
appropriate when it helps resolve the design. Once the design is accepted, straightforward fixes within
that design proceed directly to implementation. For routine changes within agreed architecture, implement
the smallest coherent change without a separate design phase.

Format changed C++ before review. Validate changed behavior with the smallest meaningful set of checks:
reuse or adapt existing tests before adding cases, and add coverage only for a distinct failure risk or
contract not already protected. Update the owning user documentation when its contract changes. Remove a
TODO only when its complete deliverable is qualified. Stop for user approval before out-of-scope public
API changes, automatic chemistry policy, heavyweight dependencies, destructive commands, or broad rewrites.
Preserve the scientific and architectural rules above.

### Test design

- Test observable behavior, not a second copy of the implementation. Do not mirror production catalogs,
  metadata tables, export lists, or private file layouts in expected-value inventories. Prefer a few
  contrasting examples and behavioral invariants; adding a catalog entry should not require updating a
  parallel test catalog.
- Keep independent analytical or published numerical checks in the native scientific layer. Interface
  tests should protect their own risks, such as conversion, ownership, mapping, installation, and failure
  handling, rather than repeat the numerical suite through every entry point. Test count and coverage
  percentage are not reasons to add cases.
- Treat machine-readable schemas and documented output semantics as contracts. Check schema conformance,
  not how the schema is written. For the CLI, check exit statuses, file effects, and structured output;
  do not parse human-facing reports, help, progress rendering, or diagnostic prose. Avoid exact `repr`
  snapshots. For exceptions, prefer type and structured details; check essential diagnostic context only
  when it is part of the behavior under test, not complete wording.
- Use assertions, matchers, and fixture facilities already provided by the test framework when they
  simplify the code; verify availability in the version in use. Do not reimplement exception-checking or
  failure-reporting machinery that the framework already supplies.
- Share genuinely identical, substantial setup in existing test support, keeping scenario inputs and
  assertions visible. Prefer a local case loop for repeated scenarios that differ only in input values.
  Keep deliberately different scientific fixtures separate. Avoid pass-through wrappers, configurable
  test frameworks, and fixture classes that add more indirection than duplication they remove. A test
  cleanup should reduce maintenance burden, not grow the suite to compensate for removed snapshots.

### Validation cadence

- After design is stable, build the affected `gcc-debug` targets and run their focused tests.
  Review the architectural diff before expensive qualification. Run the full `gcc-debug` test suite once
  for a coherent feature that is ready to commit, not for every edit or review correction. Documentation-
  only changes need relevant content/link checks, not builds; check executable examples when affected.
- Select additional profiles by risk: `clang-debug` for compiler-sensitive behavior such as templates,
  conversions, or overloads; release profiles for numerical methods, Eigen, optimization-, or
  `NDEBUG`-sensitive behavior; `clang-asan` for lifetime, bounds, mapping, parser, container, view, or
  pointer risks; and `clang-ubsan` for arithmetic, conversion, alignment, indexing, or other undefined-
  behavior risks. Run focused checks after design stabilizes; run sanitizer workflows sequentially with
  conservative parallelism.
- Run `clang-tidy` on affected implementation at qualification, not during sketch iterations; tests and
  profiles establish delivered behavior, so speculative drafts do not need full qualification. A broad
  debug, release, sanitizer, and clang-tidy matrix is for major integration or release qualification, or
  when requested; applicable CI evidence may supply it. Fixes rerun affected checks and do not
  automatically trigger the full matrix. Report failures rather than hiding them.

The executable commands for every preset and check are maintained in
[DEVELOPMENT.md](DEVELOPMENT.md).

Do not run the CLI directly from the build tree for parameter-dependent commands: bundled parameter
sets are installed resources. Install locally and run the installed executable, or use the registered
CTest CLI tests, which prepare an installation.

### Container validation

Use the distributable CLI image and compatibility-container workflows in
[DEVELOPMENT.md](DEVELOPMENT.md#installation-and-containers). Run compatibility builds one at a time and
update pinned distribution releases intentionally.

## C++ conventions

- Use C++23. Public names remain under
  `chargefw::{core,features,parameters,methods,charges,calculation,adapters}`.
- Follow `.clang-format`: LLVM-derived, four spaces, 100-column limit, left-attached pointers and
  references.
- Prefer value ownership, `std::span` for read-only views, and explicit non-owning references/pointers
  when caller lifetime is guaranteed.
- Use exceptions for invalid inputs and calculation failures; preserve actionable target context.
