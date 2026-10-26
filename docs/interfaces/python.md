# Python

The `chargefw` package reads molecules, calculates charges, and returns them as NumPy arrays.
[Installation](../getting-started/installation.md#python-package) describes how to install it, and the
[examples](../examples/index.md) show complete scripts.

```python
import chargefw

molecules = chargefw.io.read("molecules.sdf", format="sdf")
result = chargefw.calculate(molecules, method="eem", parameter_set="EEM_Cheminf_b3lyp_npa")
print(result.assignments[0].values)
```

Calculations run in the calling thread, but release Python's global interpreter lock while they work,
so several calculations can run at the same time in different threads.

## Reading molecules

`chargefw.io.read()` reads a file and `chargefw.io.parse()` reads text. Both need the format:
`"mol"`, `"sdf"`, `"mol2"`, `"molecule-json"`, `"pdb"`, or `"mmcif"`; the file extension is not used.
`chargefw.io.INPUT_FORMATS` lists them. Both read the whole input and return a `MoleculeCollection`:

| Input | Molecules in the collection | Options |
| --- | --- | --- |
| MOL | One molecule with one conformer | |
| SDF, MOL2 | One molecule per record, one conformer each | |
| Molecule JSON | One molecule per entry | `conformers` |
| PDB | One molecule; models that match the first become conformers | `selection`, `bonds`, `conformers` |
| mmCIF | One molecule per data block with coordinates | `selection`, `bonds`, `conformers` |

For PDB and mmCIF, `selection` is `"all"` (default), `"polymers-and-ligands"`, or `"polymers"`;
`bonds` is `"hybrid"` (default), `"templates"`, `"explicit"`, or `"none"`; `conformers` is `"all"`
(default) or `"first"`. The [format reference](../reference/formats.md#pdb-and-mmcif-input) explains them,
and [Preparing input](../guide/input.md) what a molecule needs before a calculation.

```python
molecules = chargefw.io.read("structure.cif", format="mmcif", selection="polymers-and-ligands")

for molecule in molecules:
    print(molecule.name, molecule.atom_count, molecule.bond_count,
          molecule.conformer_names, int(molecule.formal_charges.sum()))
```

### Molecules from arrays

Build a `Molecule` directly when you have the data as arrays:

```python
molecule = chargefw.Molecule(
    atomic_numbers=[8, 1, 1],
    formal_charges=[0, 0, 0],
    bonds=[[0, 1, 1], [0, 2, 1]],
    coordinates=[[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]],
    name="water",
)
```

| Argument | Value |
| --- | --- |
| `atomic_numbers` | Integers of shape `(N,)`, 1 to 100 |
| `formal_charges` | Integers of shape `(N,)`; default all zero |
| `bonds` | Rows `(first_atom, second_atom, order)` of shape `(B, 3)`, with order 1, 2, or 3 |
| `coordinates` | Ångströms, shape `(N, 3)` for one conformer or `(C, N, 3)` for several |
| `name`, `atom_names`, `conformer_names` | Optional names |
| `atom_ids` | Optional IDs, strings or 64-bit integers; default the atom positions |
| `source`, or `source_name`, `record_index`, `record_id` | Optional record of where the molecule came from |

ChargeFW copies the arrays and checks them: bonds of an atom to itself, repeated bonds, indices out of
range, non-integer values, and shapes that do not fit are errors. The molecule's arrays are read-only,
and `coordinates` always has the shape `(C, N, 3)`, with `C = 0` when there are none. Coordinates that
are not finite are kept, but methods that need geometry then reject the molecule.

`chargefw.MoleculeCollection` holds several molecules, in order, with an optional name. The molecules
and collections cannot be changed after they are created.

### Where atoms came from

Molecules read from a file have a `source_mapping` that links them back to the file:

- `source_mapping.atoms`: for every atom, its position in the file and its ID where the format has
  IDs;
- `source_mapping.conformers`: for every conformer, its model or record and the file sites of its atoms,
  with chain, residue, and atom labels, alternate location, and insertion code for PDB and mmCIF;
- `source_mapping.components`: for PDB and mmCIF, the residues, each with its name and atom indices;
- the format and the reading options that were used, and whether the file contained bonds.

A molecule built from arrays has `source_mapping` set to `None`, even when the arrays came from another
molecule. `molecule.atom_ids` holds the file's atom IDs, or the atom positions when the format has none.

### Periodic table

`chargefw.periodic_table` holds the elements by atomic number. Index it by number or by symbol to get an
`Element` with `atomic_number`, `symbol`, `name`, `period`, and `group`; `group` is 1 to 18, or `None`
for lanthanides and actinides:

```python
symbols = [chargefw.periodic_table[number].symbol for number in molecule.atomic_numbers]
print(chargefw.periodic_table["Se"].atomic_number, "Fm" in chargefw.periodic_table)
```

## Calculating charges

`chargefw.calculate()` takes a molecule, a collection, or any iterable of molecules, and the following
keyword arguments:

| Argument | Value |
| --- | --- |
| `method` | Method ID or `Method`; default `None`, chosen automatically |
| `parameter_set` | Parameter-set ID or `ParameterSet`; requires `method` |
| `options` | Options of the named method, such as `{"iters": 8}`; requires `method` |
| `options_by_method` | Options per method ID for automatic selection, such as `{"peoe": {"iters": 8}}` |
| `parameter_matching` | `"strict"` (default) or `"permissive"`; see [permissive matching](../reference/parameter-format.md#strict-and-permissive-matching) |
| `execution` | `"auto"` (default), `"full"`, `"cutoff"`, or `"cover"`; see [Large structures](../guide/large-structures.md) |
| `radius` | Fragment radius in ångströms, at least `8.0`; required with `"cutoff"` and `"cover"`, optional with `"auto"` (default 12), not allowed with `"full"` |
| `cutoff_threshold` | Size above which `"auto"` uses cutoff; default `20_000` atoms, `None` for no limit |
| `cover_threshold` | Size above which `"auto"` uses cover; default `80_000` atoms, `None` for no limit |
| `threads` | Maximum number of threads; by default, and with `0`, all available cores |
| `fixed_ions` | Residue names of ions to keep as fixed charges; see [below](#fixed-ions) |
| `observer` | Receives progress and can cancel; see [Progress and cancellation](#progress-and-cancellation) |

`options` and `options_by_method` cannot be used together. Thresholds must satisfy
`cutoff_threshold <= cover_threshold`. To remove both limits, pass `cutoff_threshold=None` and
`cover_threshold=None`; `cutoff_threshold=None` alone is an error, `cover_threshold=None` alone keeps the
cutoff limit.

Name both the method and the parameter set when results must be reproducible. Automatic selection, as
described in [Choosing a method](../guide/methods.md#choosing-a-method), is deterministic, but its order
is not a recommendation for a particular molecule.

### One model for all molecules

A single `calculate()` call uses one method, parameter set, and execution mode for all molecules it is
given. If they do not apply to one molecule of a collection, the calculation fails. Automatic selection
likewise chooses a model that applies to every molecule. `result.assignments_by_molecule` groups the
charges by input molecule:

```python
result = chargefw.calculate(molecules, method="eem", parameter_set="EEM_Cheminf_b3lyp_npa")

for molecule, assignments in zip(molecules, result.assignments_by_molecule, strict=True):
    for assignment in assignments:
        print(molecule.name, assignment.conformer_index, assignment.values)
```

To let every molecule get its own model, calculate them one by one; each call then selects
separately:

```python
for molecule in molecules:
    result = chargefw.calculate(molecule)
    print(molecule.name, result.plan.method.id, result.plan.policy.mode)
```

### Fixed ions

`fixed_ions` takes residue names such as `["MG", "CA"]`, or `chargefw.COMMON_IONS` or
`chargefw.ALL_IONS`. Each selected residue must be a single atom of the ion's element, without bonds;
names that are not in the structure are ignored. [Fixed ions](../guide/fixed-ions.md) describes the
model and the catalog.

```python
result = chargefw.calculate(molecules, method="sqeqp", parameter_set="SQEqp_Schindler2021_CCD_gen",
                            fixed_ions=["MG"])
```

## Results

`calculate()` returns a `CalculationResult`:

| Attribute | Content |
| --- | --- |
| `assignments` | All sets of charges, as `ChargeAssignment` values |
| `assignments_by_molecule` | The same, grouped by input molecule |
| `plan` | The method, parameter set, options, execution mode and radius, and warnings that were used; a record, not a plan you can calculate again (use [`assess()`](#assessing-before-calculating) for that) |
| `requested` | What was requested |
| `molecules` | The input molecules |
| `status` | `"success"`, or the reason for a failure |
| `rejections`, `warnings` | Rejected methods and parameter sets with reasons, and warnings |
| `failure_message` | Text of a failure, if any |
| `timings` | Time spent selecting and calculating |

A `ChargeAssignment` has the charges in `values`, a read-only NumPy `float64` array in atom order, and
`molecule_index`, `conformer_index` (`None` for methods without geometry), `source`, and `atom_ids`.
`numpy.asarray(assignment)` also returns the charges. [Results](../guide/results.md) explains
assignments, conformers, and atom order.

`result.assignment(molecule=0, conformer=1)` returns one assignment. For methods that use geometry,
`conformer` is required, also for molecules with one conformer; for the other methods, omit it.

With fixed ions, `result.requested.fixed_ions` holds the names you passed and `result.plan.fixed_ions`
the ions that were fixed, each with its molecule and atom index and its charge.

### Errors

Invalid arguments raise `TypeError` or `ValueError`. A calculation that fails raises a
`chargefw.ChargeFWError`:

- `NoExecutablePlanError`: no method and parameter set apply with the requested settings;
- `NumericalFailureError`: the calculation failed, for example did not converge;
- `CalculationCancelledError`: the calculation was cancelled.

The exception's `result` holds the complete result, with the reasons; it can still be written as result
JSON:

```python
try:
    result = chargefw.calculate(molecules, method="qeq", parameter_set="QEq_original")
except chargefw.ChargeFWError as error:
    chargefw.io.write("failed-result.json", error.result, format="result-json")
    raise
```

Indices of molecules, atoms, bonds, and conformers in the reasons count from zero.

## Writing results

`chargefw.io.write()` writes a result to a file and `chargefw.io.dumps()` returns it as text. The format
is `"result-json"`, `"mol2"`, or `"mmcif"`; `chargefw.io.OUTPUT_FORMATS` lists them:

```python
chargefw.io.write("charged.cif", result, format="mmcif")
chargefw.io.write("charged.mol2", result, format="mol2")
```

MOL2 and mmCIF can be written for input in any format. [Choosing an output
format](../guide/results.md#choosing-an-output-format) compares the formats, and the
[format reference](../reference/formats.md#charge-output) describes them. An unsupported format or a
result the format cannot hold raises `ValueError`, and a file that cannot be written raises an `OSError`,
such as `PermissionError`. A failed write can leave an incomplete file, but the result itself stays
usable.

Record and atom IDs are strings or 64-bit integers; NumPy integers are converted to Python integers.
Result JSON keeps integer IDs as numbers. MOL2 and mmCIF write them as text where the format needs a name.

## Methods and parameter sets

`chargefw.methods` and `chargefw.parameter_sets` look up methods and parameter sets by ID. Methods are
ordered by descending priority, then ID, as in the [method table](../guide/methods.md#available-methods);
parameter sets are grouped by method ID and, within a method, ordered by descending priority, then ID.
Pass the ID or the object to `calculate()`:

```python
eem = chargefw.methods["eem"]
print(eem.full_name, eem.time_complexity, eem.requires_coordinates, eem.supports_cover)

for parameter_set in eem.parameter_sets.values():
    print(parameter_set.id, parameter_set.name)

iterations = chargefw.methods["peoe"].options["iters"]
print(iterations.default, iterations.minimum, iterations.maximum)
```

A `Method` has its name, publication, `notes`, `priority`, `requires_coordinates`,
`time_complexity` and `memory_complexity`, `supports_cutoff`, `supports_cover`,
`supports_fixed_point_sources`, `options`, and `parameter_sets`. A method that requires coordinates
needs finite coordinates for every conformer and gives one assignment per conformer. These attributes
describe the method in general; whether it applies to particular molecules shows only when you
calculate or [assess](#assessing-before-calculating). `formal` and `dummy` are not charge models: they
return the formal charges and zeros.

A `ParameterSet` has its metadata and `atom_keys` and `bond_keys`, the atom and bond types it has values
for, in the order they are matched:

```python
peptide = chargefw.parameter_sets["SQEqp_Schindler2021_PUB_pept"]
for key in peptide.atom_keys:
    print(key.element, key.classifier, key.type)
for key in peptide.bond_keys:
    print([(atom.element, atom.type) for atom in key.atoms], key.classifier, key.type)
```

Element `"*"` matches every element. The [parameter-set format](../reference/parameter-format.md)
explains classifiers and types. Python offers the parameter sets installed with ChargeFW; it cannot load
other parameter files.

## Assessing before calculating

`chargefw.assess()` takes the same molecules and arguments as `calculate()`, without `observer`, and
checks which methods and parameter sets apply, without calculating. Use it to see the reasons for
rejections, compare the possible plans, or calculate several plans without preparing the molecules
again; the [parameter-set comparison example](../examples/compare_parameter_sets.md) does that.

```python
assessment = chargefw.assess(molecule, method="eem", execution="full")

for plan in assessment.plans:
    print(plan.method.id, plan.parameter_set.id, plan.policy.mode, plan.warnings)

if assessment.default_plan is None:
    for rejection in assessment.rejections:
        for issue in rejection.issues:
            print(rejection.method.id, issue.kind, issue.message)
```

An `Assessment` has the possible `plans` in selection order, the `rejections`, the `default_plan` that
`calculate()` would use, and the time taken in `seconds`. A `Plan` has its method, parameter set, all
option values, execution mode, and warnings. `plan.parameter_set` is `None` for methods without
parameter sets.

Calculate a plan by passing it after the molecules; the molecules are not prepared again. A plan
estimates its memory and number of operations before you run it:

```python
plan = assessment.default_plan
if plan is not None:
    estimate = plan.estimate(threads=4)
    print(estimate.memory_bytes / 2**30, estimate.operations)
    result = chargefw.calculate(molecule, plan, threads=4)
```

[Threads and estimates](../guide/large-structures.md#threads-and-estimates) explains the numbers.
`threads=None` in `estimate()` uses the `threads` of the assessment, and `threads` in `calculate()`
overrides it for one calculation.

A plan works only with the molecules it was assessed with: the same molecule objects, in the same
order, in a collection with the same name; copies or other molecules are rejected. Selection arguments,
including `fixed_ions`, cannot be given together with a plan. A plan
stays usable after the assessment is deleted, and several threads can calculate the same plan at once.

## Progress and cancellation

Pass a subclass of `chargefw.CalculationObserver` as `observer` to receive progress and to cancel:

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
result = chargefw.calculate(molecule, method="qeq", parameter_set="QEq_original", observer=observer)
```

Each `CalculationProgress` has a `phase`: `"computation_started"`, `"target_started"`,
`"fragment_progress"`, `"target_finished"`, or `"computation_finished"`. A target is one conformer of a
molecule for methods that use geometry, and one molecule for the others. Fragment progress comes only
from cutoff and cover, at limited frequency. The event also has the execution mode, method ID, target
and fragment counts, molecule and conformer index, and elapsed seconds; fields unrelated to the phase are
zero or `None`.

`on_progress()` can be called from several worker threads, in any order; it holds the global interpreter
lock, but protect shared state with a lock as above. An exception in it is reported through
`sys.unraisablehook` and does not stop the calculation. Selection before the calculation sends no
progress.

When `cancelled()` returns true, ChargeFW stops before the next molecule, conformer, or fragment and
raises `CalculationCancelledError`; its result has status `"cancelled"` and no charges.

## RDKit

`pip install --pre "chargefw[rdkit]"` installs the RDKit integration in `chargefw.io.rdkit`.

```python
from chargefw.io import rdkit as chargefw_rdkit

molecule = chargefw_rdkit.from_mol(rdkit_molecule, bond_conversion="single")
result = chargefw.calculate(molecule, method="qeq", parameter_set="QEq_original")
chargefw_rdkit.attach_charges(rdkit_molecule, result)
```

`from_mol()` copies an RDKit molecule as it is: it does not sanitize, add or remove hydrogens, protonate,
or generate coordinates. Every hydrogen must be an explicit atom; call `Chem.AddHs()` first for molecules
from SMILES. By default, only single, double, and triple bonds are accepted. `bond_conversion="single"`
turns aromatic and dative bonds into single bonds, as ChargeFW's file readers do with aromatic bonds;
other bond types, such as ionic, hydrogen, or zero-order bonds, remain errors.

`attach_charges()` stores the charges of one assignment as a double property of each atom, named by
`property_name` (default `ChargeFWPartialCharge`), and as a property list that RDKit writes to SD files.
Existing properties of that name are an error unless you pass `overwrite=True`. It matches atoms by their IDs, which must
be integers that map every atom exactly once. `molecule_index` selects the input molecule (default `0`).
When that molecule has a single assignment, it is used; otherwise pass `conformer`, the ChargeFW
conformer index, not the RDKit conformer ID. The property belongs to the molecule, not to a conformer, so
keep separate copies or property names to store several sets of charges.

## Gemmi

ChargeFW reads PDB and mmCIF files without the Gemmi Python package. To convert Gemmi objects, install
the integration with `pip install --pre "chargefw[gemmi]"`. `chargefw.io.gemmi` can be imported
without it; its functions then raise an `ImportError` that says what to install.

```python
import chargefw.io.gemmi

molecules = chargefw.io.gemmi.from_document(document)  # or from_structure(structure)
result = chargefw.calculate(molecules, method="qeq", parameter_set="QEq_original")

charged = chargefw.io.gemmi.to_document(result)
charged.write_file("charged.cif")

# Or add the charges to the original document, keeping all its other data.
chargefw.io.gemmi.attach_charges(document, result)
document.write_file("charged-with-source-metadata.cif")
```

`from_structure()` and `from_document()` write the Gemmi object as mmCIF text and read it with
ChargeFW's own reader, so the `selection`, `bonds`, and `conformers` options behave as for files.

`to_document()` creates a new mmCIF document, the same as `chargefw.io.dumps(..., format="mmcif")`; it
does not copy other data of the input.

`attach_charges(document, result)` adds the charges to the document that was passed to
`from_document()`. It first checks that the document still matches the calculated molecules: blocks,
atom sites, IDs, models, elements, labels, coordinates, and formal charges. Occupancies and B factors may
differ. If anything does not match, the document stays unchanged. Existing charge data is an error unless
you pass `overwrite=True`. The function replaces the document's blocks, so get block objects such as
`document.sole_block()` again afterwards.
