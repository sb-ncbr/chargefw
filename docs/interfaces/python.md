# Python package

The `chargefw` Python package provides owned molecule values, immutable NumPy arrays, bundled method and
parameter catalogs, direct calculation, optional reusable assessment plans, calculation results, and
Gemmi conversion over the native ChargeFW engine.

The package is synchronous and in-process. Native molecule construction, assessment, and calculation
release the GIL, and independent calculations can run concurrently. Calculations can report structured
progress and support cooperative cancellation through a per-request observer.

The [examples](../examples/index.md) show complete workflows, including toolkit integration.

See [Installation](../getting-started/installation.md#python-package) for installing the package.

## Quick start

```python
import chargefw

molecules = chargefw.io.read("structure.sdf", format="sdf")
result = chargefw.calculate(
    molecules,
    method="qeq",
    parameter_set="QEq_original",
    execution="full",
)
print(result.assignments[0].values)
```

For reproducible scientific work, select both the method and its parameter set when the method uses one.
Automatic selection is deterministic, but catalog priority is not a recommendation that one charge model
or parameterization is scientifically preferable for a particular molecule.

## Reading molecular files

Serialized molecular formats are read through `chargefw.io`. `parse()` accepts text, while `read()`
accepts string or path-like filesystem paths:

```python
import chargefw

molecules = chargefw.io.read(
    "structure.sdf",
    format="sdf",
)
```

`parse()` and `read()` require an explicit `format` selected from `"mol"`, `"sdf"`, `"mol2"`,
`"molecule-json"`, `"pdb"`, and `"mmcif"`. `parse()` also accepts an optional `source_name`. Both return
a `MoleculeCollection`, including formats that contain exactly one molecule. File extensions are not
inspected, and the current Python API reads and materializes the complete collection eagerly. The
supported values are available as `chargefw.io.INPUT_FORMATS`.

The Python-facing collection distinctions are:

| Input | Python collection behavior | Python-specific choices |
| --- | --- | --- |
| MOL | One molecule with one conformer | None |
| SDF and MOL2 | One molecule per record, one conformer each | None |
| Molecule JSON | One molecule per array entry | `conformers="first"` or `"all"` |
| PDB | One molecule; compatible models become conformers | `selection`, `bonds`, `conformers` |
| mmCIF | One molecule per coordinate-bearing block | `selection`, `bonds`, `conformers` |

Exact parsing, normalization, diagnostics, and record semantics belong to the
[molecular format reference](../reference/formats.md#format-overview).

Inspect imported molecules through immutable arrays and source identities:

```python
for molecule in molecules:
    print(
        molecule.source,
        molecule.name,
        molecule.atom_count,
        molecule.bond_count,
        molecule.conformer_names,
        int(molecule.formal_charges.sum()),
    )
```

`calculate()` accepts a `MoleculeCollection` or any iterable of `Molecule` values. One calculation plan
applies to the complete collection. If a method or parameter set is inapplicable to one SDF record or
mmCIF block, it is not executable for that collection. Automatic selection also chooses one method,
parameter set, and execution mode that applies to every record.

For a collection, `assignments_by_molecule` groups the source-ordered assignments for each input molecule.
Use it with the original collection for ordinary molecule-by-molecule iteration:

```python
molecules = chargefw.MoleculeCollection([first_molecule, second_molecule])
result = chargefw.calculate(
    molecules,
    method="qeq",
    parameter_set="QEq_original",
    execution="full",
)

for molecule, assignments in zip(molecules, result.assignments_by_molecule, strict=True):
    for assignment in assignments:
        print(molecule.name, assignment.conformer_index, assignment.values)
```

Process records separately when independent per-record policy is intentional. Each call below performs
its own automatic assessment, so records can select different methods, parameter sets, or execution
modes:

```python
for molecule in molecules:
    result = chargefw.calculate(molecule)
    parameter_set = result.plan.parameter_set
    print(
        molecule.name,
        result.plan.method.id,
        parameter_set.id if parameter_set is not None else None,
        result.plan.policy.mode,
    )
```

This loop is not equivalent to `calculate(molecules)`: the latter requires one common plan and produces
one result containing source-aligned assignments for the complete collection. Choose the collection
boundary deliberately when using automatic selection.

Retrieve one assignment directly by its zero-based source indices. For geometry-dependent methods,
`conformer` is required even when the molecule has only one conformer; use `conformer=0` in that case.
For geometry-independent methods, omit it: `result.assignment(molecule=0)` returns the molecule's
single assignment.

```python
second_conformer_charges = result.assignment(molecule=0, conformer=1).values
```

## Constructing Molecule objects from arrays

Use `Molecule` directly when molecular data is already available as arrays:

```python
molecule = chargefw.Molecule(
    atomic_numbers=[8, 1, 1],
    formal_charges=[0, 0, 0],
    bonds=[[0, 1, 1], [0, 2, 1]],
    coordinates=[
        [0.0, 0.0, 0.0],
        [0.96, 0.0, 0.0],
        [-0.24, 0.93, 0.0],
    ],
    name="water",
)
```

`Molecule` accepts:

- `atomic_numbers`: an integer vector of shape `(N,)`;
- `formal_charges`: an optional integer vector of shape `(N,)`, defaulting to zero;
- `bonds`: optional integer rows `(first_atom, second_atom, order)` of shape `(B, 3)`;
- `coordinates`: `(N, 3)` for one conformer or `(C, N, 3)` for multiple conformers;
- optional molecule, atom, and conformer names;
- a `SourceIdentity` or its `source_name`, `record_index`, and `record_id` fields; and
- optional string or signed 64-bit integer `atom_ids` for round-trip mapping.

Inputs are validated, normalized to C-contiguous `int64` or `float64` arrays, and copied. Public arrays
are read-only, so later mutation or destruction of the caller's arrays cannot change the molecule.
Coordinates are always exposed with shape `(C, N, 3)`; `None` and `(0, N, 3)` both mean no conformers.
Coordinates are interpreted in ångströms; the Python API does not attach units or convert them.

Molecules returned by serialized readers expose an immutable `source_mapping`. Its `atoms` tuple maps
calculation atom order to zero-based source positions and exact source ID tokens where the format provides
them. Its `conformers` tuple retains each selected source conformer and its calculation-ordered source
sites. Structural sites additionally expose model identity, selected alternate location, entity, insertion
code, segment where available, and separate mmCIF author and label atom/residue/chain/sequence identities.
The mapping also reports source format, record-local selection/conformer/bond policy, alternate-location
policy, and whether source connectivity was absent, explicitly empty, or present. Imported `atom_ids` use
exact source IDs when present and source positions otherwise. A newly constructed molecule, including one
built from an imported molecule's arrays, has `source_mapping is None` and does not inherit verified
correspondence. Structural imports expose their selected residue partition as immutable
`source_mapping.components` entries (`SourceComponentInstance`), each with a canonical `component_id` and
calculation-ordered `atom_indices`; see the [format reference](../reference/formats.md#pdb-and-mmcif-input) for its
ordering and conformer contract.

Atomic numbers 1–100 are accepted because they are represented by the bundled periodic table; individual
methods and parameter sets may support a smaller subset. Only bond orders 1, 2, and 3 are accepted. Self
bonds, duplicate bonds, invalid indices, non-integral integer input, unsupported atomic numbers, and
mismatched shapes are rejected. Non-finite coordinates are retained, but methods that require geometry
are inapplicable when coordinates are missing or non-finite.

`chargefw.periodic_table` is the immutable bundled table in atomic-number order. Index it by atomic
number, including NumPy integers, or by exact element symbol to get an `Element` with `atomic_number`,
`symbol`, `name`, `period`, and `group`. `group` is the IUPAC group 1–18, and `None` for lanthanides and
actinides placed in the f-block:

```python
symbols = [chargefw.periodic_table[number].symbol for number in molecule.atomic_numbers]
print(chargefw.periodic_table["Se"].atomic_number, "Fm" in chargefw.periodic_table)
```

`MoleculeCollection` is an immutable source-ordered sequence with an optional name. Omitted atom IDs
default to source-order integers.

## Methods and parameter sets

The package-level `chargefw.methods` and `chargefw.parameter_sets` values are immutable ordered mappings
with lookup by stable ID. Select by ID or pass a descriptor from these mappings to `assess()` or
`calculate()`:

```python
eem = chargefw.methods["eem"]

for method_id, method in chargefw.methods.items():
    print(method_id, method.full_name)

print(eem.time_complexity, eem.memory_complexity)
print(eem.requires_coordinates, eem.supports_cutoff, eem.supports_cover)

for parameter_set in eem.parameter_sets.values():
    print(parameter_set.id, parameter_set.name)

iterations = chargefw.methods["peoe"].options["iters"]
print(iterations.default, iterations.minimum, iterations.maximum)

peptide = chargefw.parameter_sets["SQEqp_Schindler2021_PUB_pept"]
for key in peptide.atom_keys:
    print(key.element, key.classifier, key.type)
for key in peptide.bond_keys:
    print([(atom.element, atom.type) for atom in key.atoms], key.classifier, key.type)
```

Method descriptors expose names, publication metadata, human-readable `notes`, priority, coordinate
requirements, Big-O time and memory complexity, reduced-mode capabilities, fixed point source support,
options, and associated parameter sets. `supports_fixed_point_sources` reports whether the method can
numerically account for fixed external point charges; it does not establish applicability for a particular
molecule and parameter set. Notes describe important implementation scope or input behavior; assessment
remains the authoritative applicability check for a specific request.
`requires_coordinates=True` means every input molecule must have coordinates and every conformer must
contain finite, non-coincident coordinates; the method then produces one assignment per conformer.
Complexity uses `n` for atoms and `m` for bonds, as defined in the
[Methods and parameter sets](../guide/methods.md). The `formal` method copies input formal
charges and `dummy` returns zeros; neither is an empirical partial-charge model.

Parameter-set descriptors expose metadata and their entry keys, which show the atom and bond
environments a set parameterizes separately. `atom_keys` is a tuple of
`AtomParameterKey(element, classifier, type)` values; element `"*"` matches every element.
`bond_keys` is a tuple of `BondParameterKey(atoms, classifier, type)` values, where `atoms` holds the two
endpoint atom keys and the classifier and type describe the bond. Both tuples follow the file's entry
order, which is the matching precedence. Python exposes only the installed bundled parameter catalog.
Its package-level and method-level mappings are immutable and cannot be extended. The
[parameter-set JSON reference](../reference/parameter-format.md) documents the bundled data and classifier behavior,
including permissive matching. Selection order and method scope are defined in
[Methods and parameter sets](../guide/methods.md).

## Calculation policy

Direct `calculate()` accepts keyword-only policy arguments:

| Argument | Values and defaults |
| --- | --- |
| `method` | Method ID, `Method`, or `None` |
| `parameter_set` | Parameter-set ID or `ParameterSet`; requires `method` |
| `options` | Option mapping for an explicitly selected method |
| `options_by_method` | Method-ID to option-mapping overrides for automatic selection |
| `parameter_matching` | `"strict"` (default) or `"permissive"` |
| `execution` | `"auto"` (default), `"full"`, `"cutoff"`, or `"cover"` |
| `radius` | Reduced radius in ångströms; explicit cutoff/cover require at least `8.0` |
| `cutoff_threshold` | Automatic full-to-cutoff threshold; default `20_000`, `None` is unlimited |
| `cover_threshold` | Automatic cutoff-to-cover threshold; default `80_000`, `None` is unlimited |
| `threads` | Non-negative oneTBB thread limit; omitted or `0` delegates to oneTBB |
| `fixed_ions` | Sequence of named imported monatomic ion components with no graph bonds; default `None` |

Finite thresholds are non-negative integers with `cutoff_threshold <= cover_threshold`. To make both
thresholds unlimited, pass `cutoff_threshold=None, cover_threshold=None`; setting only
`cutoff_threshold=None` leaves the default finite cover threshold and raises `ValueError`. Setting only
`cover_threshold=None` is valid and keeps the default cutoff threshold.

The [Large structures](../guide/large-structures.md) defines the shared execution modes, automatic
selection policy, conservation behavior, and approximation limits.

Pass the same named ion selection to `assess()` or direct `calculate()` to apply fixed-ion handling. Each
selected component must contain one atom of its catalog element and have no graph bonds. The request
snapshots names to an immutable tuple; `None` and an empty sequence disable fixed-ion handling.
The [Fixed ions](../guide/fixed-ions.md) describes how imported components are resolved,
which methods support the model, and its active-charge contract.

`chargefw.COMMON_IONS` and `chargefw.ALL_IONS` are immutable tuples of component names accepted by the
same argument. `COMMON_IONS` selects the common-ion subset; `ALL_IONS` selects every monatomic ion in the
[bundled catalog](../guide/fixed-ions.md). Both expand explicit caller intent, using bundled
template charges; only matching components present in the input become fixed sources.

```python
assessment = chargefw.assess(
    molecules,
    method="sqeqp",
    parameter_set="SQEqp_Schindler2021_CCD_gen",
    execution="full",
    fixed_ions=["MG"],
)
result = chargefw.calculate(molecules, assessment.default_plan)

common = chargefw.calculate(molecules, fixed_ions=chargefw.COMMON_IONS)
all_ions = chargefw.calculate(molecules, fixed_ions=chargefw.ALL_IONS)
```

Flat `options` require an explicit method. `options` and `options_by_method` cannot be combined. Automatic
execution accepts an optional radius override. Explicit full execution rejects a radius; explicit cutoff
and cover require one. Direct calculation performs assessment internally, executes the selected default
plan, and returns its effective policy and warnings in the result.

```python
result = chargefw.calculate(
    molecule,
    method="eem",
    parameter_set="EEM_Baek1991",
    execution="auto",
)
print(result.plan.method.id, result.plan.policy.mode, result.warnings)
```

## Advanced assessment and plan reuse

Calling `assess()` directly is optional. Use it to inspect applicability or policy warnings before
calculation, compare runnable alternatives, or execute several plans without repeating molecule
preparation and parameter classification. The
[parameter-set comparison example](../examples/compare_parameter_sets.md) demonstrates this workflow with one
assessment and several plans.

`assess()` accepts the same selection and policy arguments as direct `calculate()`. Its `threads`
argument configures default plan executions; passing `threads` to `calculate(..., plan, threads=...)`
overrides that default for one execution.

```python
assessment = chargefw.assess(molecule, method="eem", execution="full")

for plan in assessment.plans:
    parameter_set = plan.parameter_set
    print(
        plan.method.id,
        parameter_set.id if parameter_set is not None else None,
        plan.policy.mode,
        plan.warnings,
    )
```

An `Assessment` contains priority-ordered reusable `plans`, structured `rejections`, a `default_plan`,
and assessment time in `seconds`. Each `Plan` retains the prepared native state and exposes its method,
parameter set, complete validated options, concrete execution policy, and warnings.

When no plan is executable, inspect issues without parsing human-readable exception text:

```python
if assessment.default_plan is None:
    for rejection in assessment.rejections:
        for issue in rejection.issues:
            print(rejection.method.id, issue.kind, issue.message)
```

Issue molecule, atom, bond, and conformer indices are zero-based when present.

Execute an exact assessed plan without repeating preparation or parameter classification:

```python
plan = assessment.default_plan
if plan is not None:
    result = chargefw.calculate(molecule, plan, threads=1)
```

Estimate a plan's peak working memory and operation count before executing it:

```python
estimate = plan.estimate(threads=4)
print(estimate.memory_bytes / 2**30, estimate.operations)
```

`threads=None` uses the assessment's `threads`. The returned `ResourceEstimate` is approximate; the
[resource estimates](../guide/large-structures.md#threads-and-estimates) describes what it covers.

Plans are bound to the exact molecule objects and collection name used during assessment. They cannot be
applied to another record or a reconstructed equivalent molecule, and selection arguments cannot be
supplied with a plan, including `fixed_ions`. A plan is reusable after its `Assessment`
is released and can be used by independent concurrent calculations over those same objects.

## Progress and cancellation

Subclass `CalculationObserver` and pass it to `calculate()` with the keyword-only `observer` argument:

```python
from threading import Event, Lock

import chargefw


class RecordingObserver(chargefw.CalculationObserver):
    def __init__(self) -> None:
        self.events: list[chargefw.CalculationProgress] = []
        self.cancellation = Event()
        self.lock = Lock()

    def on_progress(self, progress: chargefw.CalculationProgress) -> None:
        with self.lock:
            self.events.append(progress)

    def cancelled(self) -> bool:
        return self.cancellation.is_set()


observer = RecordingObserver()
result = chargefw.calculate(
    molecule,
    method="qeq",
    parameter_set="QEq_original",
    execution="full",
    observer=observer,
)
```

Progress phases are `"computation_started"`, `"computation_finished"`, `"target_started"`,
`"target_finished"`, and `"fragment_progress"`. A target is one molecule/conformer pair for a
geometry-dependent method and one molecule for a geometry-independent method. Fragment events are
throttled and are emitted only by cutoff and cover execution.

Each immutable `CalculationProgress` snapshot contains the phase, effective execution mode, method ID,
target indices and count, fragment completion and count, molecule and optional conformer indices, and
elapsed computation seconds. Fields unrelated to the current phase retain their default zero or `None`
value.

Callbacks may originate on oneTBB worker threads. The Python GIL is held during each callback, but
observers that update compound state should still use synchronization and must not rely on callback
arrival order during parallel execution. Callback exceptions are reported through Python's unraisable
exception hook and do not alter calculation control flow.

Returning true from `cancelled()` requests termination at the next
[cancellation check point](../guide/results.md). The call then raises
`CalculationCancelledError`; its `result` has status `"cancelled"` and contains no partial charge
assignments. Assessment itself is not observed.

## Results and failures

`CalculationResult` exposes:

- `status`;
- the immutable input collection in `molecules`;
- immutable `assignments`;
- immutable `assignments_by_molecule`, aligned with the calculation input collection;
- requested policy in `requested`;
- executed method and settings in `plan`, an `ExecutedPlan | None`;
- rejected alternatives and warnings;
- optional `failure_message`; and
- applicability and computation `timings`.

`result.plan` records the method, parameter set, options, execution policy, and warnings for the executed
plan. This `ExecutedPlan` is metadata, not a reusable `Plan`. To repeat an assessed calculation, retain a
`Plan` from `assessment.plans` or `assessment.default_plan`, as shown in
[Advanced assessment and plan reuse](#advanced-assessment-and-plan-reuse).
For fixed-ion calculations, `result.requested.fixed_ions` retains the submitted names and
`result.plan.fixed_ions` is a `FixedIons` record whose immutable `sources` tuple contains
`FixedAtomCharge` values with original molecule/atom indices and assigned charges.

Each `ChargeAssignment` contains a newly owned, read-only, C-contiguous `float64` vector together with its
molecule index, optional conformer index, `SourceIdentity`, and atom IDs.
`numpy.asarray(assignment)` returns the charge vector.

Geometry-dependent methods return one assignment per conformer. Geometry-independent methods return one
assignment per molecule without conformer identity.

Invalid Python values and request arguments raise `TypeError` or `ValueError`. Failed calculations raise
a typed `ChargeFWError` subclass:

- `NoExecutablePlanError`;
- `NumericalFailureError`;
- `CalculationCancelledError`.

Each exception retains the complete result as `exception.result`, including status, rejections, failure
text, and timings.

Failed and cancelled results can still be retained as machine-readable result JSON:

```python
try:
    result = chargefw.calculate(
        molecule,
        method="qeq",
        parameter_set="QEq_original",
        execution="full",
    )
except chargefw.ChargeFWError as error:
    chargefw.io.write("failed-result.json", error.result, format="result-json")
    raise
```

## Output

`chargefw.io.dumps()` serializes a `CalculationResult` to generated MOL2, generated mmCIF, or result JSON,
while `chargefw.io.write()` writes it as a UTF-8 file at the requested path. Both require an explicit
`format`: `"mol2"`, `"mmcif"`, or `"result-json"`. The supported values are also available as
`chargefw.io.OUTPUT_FORMATS` for discovery and argument validation. A write failure can leave the destination
incomplete.

```python
result = chargefw.calculate(
    molecules,
    method="qeq",
    parameter_set="QEq_original",
    execution="full",
)
chargefw.io.write("charged.cif", result, format="mmcif")
chargefw.io.write("charged.mol2", result, format="mol2")
```

The language-independent [charge-output reference](../reference/formats.md#charge-output) defines generated MOL2 and
mmCIF content, result JSON, coordinate requirements, precision, and mapping behavior. Python output does
not depend on the source format. Invalid argument types raise `TypeError`; unsupported output formats
and unmet format requirements raise `ValueError`. Filesystem write failures raise `OSError` subclasses,
such as `FileNotFoundError` or `PermissionError`. Output errors leave the calculation result unchanged,
so result JSON can still be serialized from the same object.

Record and atom IDs are normalized to strings or signed 64-bit integers. NumPy integer scalars are accepted
and converted to Python `int`; booleans, arbitrary hashable objects, and out-of-range integers are rejected.
Result JSON preserves integer IDs as JSON numbers rather than stringifying them. Explicit `atom_ids` from a
manually constructed molecule appear once in `input.atom_ids`. Import mapping is available through the
Python molecule and is used by output adapters.
Generated MOL2 and mmCIF use their decimal representation where a textual record or block name is
required.

## Toolkit integrations

PDB and mmCIF parsing use ChargeFW's compiled Gemmi dependency and do not require the upstream Python
package. Converting upstream `gemmi.Structure` and `gemmi.cif.Document` objects requires the optional
Gemmi Python integration:

```bash
pip install "chargefw[gemmi]"
```

Importing `chargefw.io.gemmi` remains safe without that extra; calling its conversion functions raises an
actionable `ImportError`. With the extra installed, use the object-conversion module explicitly:

```python
import chargefw.io.gemmi

structure_molecules = chargefw.io.gemmi.from_structure(structure)
document_molecules = chargefw.io.gemmi.from_document(document)
fresh_document = chargefw.io.gemmi.to_document(result)
```

The Gemmi integration serializes upstream objects through mmCIF text and then uses ChargeFW's compiled
reader. This keeps selection and bond behavior aligned with serialized input without sharing C++ objects
between extension modules.

`to_document(result)` returns a fresh `gemmi.cif.Document` through the same native writer used by
`chargefw.io.dumps(..., format="mmcif")`. It never mutates or retains an input document, and source-only
categories and original site IDs are not copied.

`attach_charges(document, result)` is the Python-only in-memory annotation path. `document` must be the
mmCIF document passed to `from_document()` for the molecules used by `result`. Before changing the live
object, ChargeFW validates coordinate-block order, block names, source site positions, exact atom IDs,
model IDs, elements, structural labels, coordinates, and formal charges against the result's owned import
mapping. Occupancy and B factors may change. Any mismatch or output error leaves the document unchanged.
Existing charge categories are rejected unless `overwrite=True`; append mode is not supported.
Successful attachment replaces the document's blocks with annotated copies, so obtain block objects such
as `document.sole_block()` again after the call. Subsequent document changes are the caller's
responsibility.

```python
result = chargefw.calculate(
    document_molecules,
    method="qeq",
    parameter_set="QEq_original",
    execution="full",
)
charged_document = chargefw.io.gemmi.to_document(result)
charged_document.write_file("charged.cif")

# Alternatively, preserve unrelated categories in the unchanged live input document.
chargefw.io.gemmi.attach_charges(document, result)
document.write_file("charged-with-source-metadata.cif")
```

Structural input defaults are `selection="all"`, `bonds="hybrid"`, and `conformers="all"`. Bond choices are
`"none"`, `"explicit"`, `"templates"`, and `"hybrid"`; selection choices are `"all"`,
`"polymers-and-ligands"`, and `"polymers"`. Their language-independent import semantics are defined in
the [PDB and mmCIF format reference](../reference/formats.md#pdb-and-mmcif-input).

RDKit is an optional dependency installed with `pip install "chargefw[rdkit]"`. The base package remains
RDKit-free. `chargefw.io.rdkit.from_mol()` copies an existing `rdkit.Chem.Mol` without sanitization,
hydrogen changes, protonation, embedding, or optimization. Every hydrogen must be an explicit RDKit atom;
call `Chem.AddHs()` first for molecules with implicit hydrogens, such as those from SMILES. By default,
it accepts only explicit RDKit single, double, and triple bond types. Pass `bond_conversion="single"` to convert aromatic bonds and the
RDKit dative bond family to native single bonds, matching the connectivity-only aromatic normalization
used by ChargeFW's serialized format readers. The option does not convert query bonds or unsupported
higher-order, fractional, ionic, hydrogen, three-center, zero-order, unspecified, or other bond types;
these remain errors. `attach_charges()` writes one selected assignment to double-valued atom properties
using integer-compatible atom IDs that map each target atom exactly once. It creates RDKit's serializable
atom-property list when that facility is available. Its `conformer` argument is a zero-based ChargeFW
conformer index, not an RDKit conformer ID, and the resulting atom property is molecule-wide rather than
conformer-scoped.

```python
from chargefw.io import rdkit as chargefw_rdkit

molecule = chargefw_rdkit.from_mol(rdkit_molecule, bond_conversion="single")
result = chargefw.calculate(
    molecule,
    method="qeq",
    parameter_set="QEq_original",
    execution="full",
)
chargefw_rdkit.attach_charges(rdkit_molecule, result)
```

The RDKit molecule in this example must already contain suitable coordinates. `attach_charges()` selects
the input molecule with `molecule_index` (default `0`). Unlike `result.assignment()`, it automatically
selects that molecule's sole assignment when `conformer` is omitted, including a single-conformer
geometry-dependent result. When the selected molecule has multiple assignments, pass `conformer=`
explicitly. Use separate molecule copies or property names if several charge vectors must be retained.
