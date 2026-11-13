# Python API reference

This page lists every public function, class, and attribute of the `chargefw` package.
[Python](../interfaces/python.md) shows how to use them together. All values are read-only, and all
indices of molecules, atoms, bonds, and conformers count from zero.

You create `Molecule`, `MoleculeCollection`, and `SourceIdentity` objects and subclasses of
`CalculationObserver`; ChargeFW creates all other objects. If you construct one of the others
yourself, for example an `AtomParameterKey` to compare with the keys of a parameter set, pass its
attributes as keyword arguments.

## Functions

| Function | Returns |
| --- | --- |
| `chargefw.calculate(molecules, *, method=None, parameter_set=None, ...)` | `CalculationResult`; the arguments are listed in [Calculating charges](../interfaces/python.md#calculating-charges) |
| `chargefw.assess(molecules, *, method=None, parameter_set=None, ...)` | `Assessment`; the same arguments as `calculate()` except `observer` |
| `chargefw.read_parameter_set(path)` | `ParameterSet` read from your own JSON file |
| `chargefw.write_parameter_set(path, parameter_set)` | Writes a parameter set as JSON |
| `chargefw.io.read(path, *, format, selection="all", bonds=None, conformers="all")` | `MoleculeCollection` read from a file |
| `chargefw.io.parse(contents, *, format, source_name="", selection="all", bonds=None, conformers="all")` | `MoleculeCollection` read from a string |
| `chargefw.io.write(path, result, *, format)` | Writes a result to a file; for mmCIF, also a list of results |
| `chargefw.io.dumps(result, *, format)` | The same output as a string |
| `chargefw.io.rdkit.from_mol(molecule, *, source_name="", bond_conversion="none")` | `Molecule` converted from an RDKit molecule |
| `chargefw.io.rdkit.attach_charges(molecule, result, *, molecule_index=0, conformer=None, property_name="ChargeFWPartialCharge", overwrite=False)` | Stores charges as RDKit atom properties |
| `chargefw.io.gemmi.from_document(document, *, source_name="", selection="all", bonds="hybrid", conformers="all")` | `MoleculeCollection` converted from a Gemmi mmCIF document |
| `chargefw.io.gemmi.from_structure(structure, *, source_name="", selection="all", bonds="hybrid", conformers="all")` | `MoleculeCollection` converted from a Gemmi structure |
| `chargefw.io.gemmi.to_document(result)` | A new Gemmi mmCIF document with the charges |
| `chargefw.io.gemmi.attach_charges(document, result, *, overwrite=False)` | Adds the charges to an existing Gemmi document |

`format` is one of `chargefw.io.INPUT_FORMATS` for reading: `"mol"`, `"sdf"`, `"mol2"`,
`"molecule-json"`, `"pdb"`, or `"mmcif"`. For writing, it is one of `chargefw.io.OUTPUT_FORMATS`:
`"mol2"`, `"mmcif"`, `"pqr"`, or `"result-json"`. `source_name` names molecules that do not come
from a file; it appears in their `source_name` and in the result JSON. `read()` uses the path.

Only mmCIF holds several results: each becomes its own charge set. The results must come from the
same molecules and use different combinations of method and parameter set; see the
[mmCIF reference](formats.md#mmcif). The other formats require exactly one result.

## Catalogs and constants

| Name | Content |
| --- | --- |
| `chargefw.methods` | Mapping of method ID to `Method`, in the order automatic selection tries them |
| `chargefw.parameter_sets` | Mapping of parameter-set ID to `ParameterSet`, grouped by method, in the order automatic selection tries them |
| `chargefw.periodic_table` | `PeriodicTable`; index it by atomic number or symbol to get an `Element` |
| `chargefw.COMMON_IONS`, `chargefw.ALL_IONS` | Residue names for `fixed_ions`; see [Fixed ions](../guide/fixed-ions.md#ion-catalog) |
| `chargefw.__version__` | Version of the package, such as `"3.0.0rc1"` |

## Molecules

### Molecule

`Molecule(atomic_numbers, formal_charges=None, bonds=None, coordinates=None, name=None, ...)` builds
a molecule from arrays; [Molecules from arrays](../interfaces/python.md#molecules-from-arrays) lists
the arguments. `len(molecule)` is its number of atoms.

| Attribute | Content |
| --- | --- |
| `atomic_numbers` | NumPy `int64` array, one value per atom |
| `formal_charges` | NumPy `int64` array, one value per atom |
| `bonds` | NumPy `int64` array of rows `(first_atom, second_atom, order)` |
| `coordinates` | NumPy `float64` array in ångströms, shape `(conformers, atoms, 3)` |
| `atom_count`, `bond_count`, `conformer_count` | Numbers of atoms, bonds, and conformers |
| `has_coordinates` | `True` when the molecule has at least one conformer |
| `name` | Name of the molecule |
| `atom_names`, `conformer_names` | Names of the atoms and conformers, such as `"CA"` or `"model-2"`; empty strings where there are none |
| `atom_ids` | ID of every atom in the file, or its position when the format has no atom IDs |
| `source_name` | File the molecule was read from, or the `source_name` you gave |
| `record_index` | Position of the record in the file, such as the third SDF record |
| `record_id` | ID of the record in the file, such as the SDF title line or the mmCIF block name; `None` when there is none |
| `source` | `SourceIdentity` with `source`, `record_index`, and `record_id` together |
| `source_mapping` | `SourceMapping` linking atoms and conformers to the file; `None` for molecules built from arrays; see [Where atoms came from](../interfaces/python.md#where-atoms-came-from) |
| `import_issues` | Tuple of `ImportIssue`: what the reader left out of the record |

### MoleculeCollection

`MoleculeCollection(molecules, name=None)` holds molecules in order. Index it and iterate over it
like a tuple; `molecules` returns the tuple.

### ImportIssue

| Attribute | Content |
| --- | --- |
| `kind` | Kind of caution, such as `"partial_charges_ignored"` or `"coordinates_2d_ignored"`; the same code appears in result JSON |
| `message` | Description |
| `line` | Line in the file, or `None` |
| `atom_indices`, `bond_indices` | Read-only NumPy `int64` arrays of the atoms and bonds the caution concerns, such as the atoms of `"unbonded_residue_atoms"` and the bonds of `"aromatic_bonds_converted"`; empty for cautions about the file |

### SourceMapping

How the atoms and conformers of a molecule correspond to its file. See
[Where atoms came from](../interfaces/python.md#where-atoms-came-from).

| Attribute | Content |
| --- | --- |
| `format` | Input format, such as `"sdf"` or `"mmcif"` |
| `atoms` | Tuple of `SourceAtomReference`, one per atom, for the first conformer |
| `conformers` | Tuple of `SourceConformerReference`, one per conformer |
| `components` | Tuple of `SourceComponentInstance`, one per residue; empty for formats without residues |
| `source_connectivity` | `"present"` when the file listed bonds, `"explicitly-empty"` when its bond section was empty, `"absent"` when it had no bond records, such as a PDB file without `CONECT`, `LINK`, or `SSBOND` records |
| `policy` | `SourcePolicy`: how the file was read |

`SourcePolicy` has the `selection`, `bonds`, and `conformers` values used for reading, as the
arguments of `chargefw.io.read()` name them, such as `"all"` and `"hybrid"`, and
`alternate_locations`, `"first-source-order"` for PDB and mmCIF, where the first
[alternate location](formats.md#structural-selection) of each atom is read. Each is `None` where the
reader has no such choice.

### SourceAtomReference

| Attribute | Content |
| --- | --- |
| `position` | Position of the atom site in the record, counting from zero |
| `id` | Atom ID in the file, such as the MOL2 atom ID or the mmCIF `_atom_site.id`; `None` for formats without atom IDs |
| `structural_labels` | `SourceStructuralLabels` for PDB and mmCIF; otherwise `None` |

### SourceConformerReference

| Attribute | Content |
| --- | --- |
| `position` | Position of the model in the file, counting from zero |
| `id` | Model number from the file, such as `"1"`; `None` when the file has none |
| `sites` | Tuple of `SourceAtomReference`: the atom sites of this model, in atom order |

### SourceStructuralLabels and SourceHierarchyLabels

| Attribute | Content |
| --- | --- |
| `author` | `SourceHierarchyLabels` with the author names: `_atom_site.auth_*` in mmCIF, the only names in PDB |
| `label` | `SourceHierarchyLabels` with the mmCIF `_atom_site.label_*` names; all `None` for PDB |
| `entity` | mmCIF entity ID, or `None` |
| `insertion_code`, `alternate_location` | As written in the file, including the mmCIF placeholders `"?"` and `"."`; `None` when absent |
| `segment` | PDB segment ID, or `None` |

`SourceHierarchyLabels` has `atom`, `residue`, `chain`, and `sequence`, each a string or `None`, such
as `"N"`, `"GLU"`, `"A"`, and `"1"`.

### SourceComponentInstance

`component_id` is the residue name, such as `"GLU"`, and `atom_indices` the indices of its atoms in the
molecule.

## Results

### CalculationResult

| Attribute | Content |
| --- | --- |
| `status` | `"success"`, `"no_executable_plan"`, `"numerical_failure"`, or `"cancelled"` |
| `assignments` | Tuple of `ChargeAssignment`; empty when the calculation failed |
| `assignments_by_molecule` | The same, as one tuple per input molecule |
| `assignment(molecule=..., conformer=None)` | One assignment; see [Selecting an assignment](#selecting-an-assignment) |
| `plan` | `ExecutedPlan` that was used, also after a numerical failure or cancellation; `None` when no plan could run |
| `requested` | `RequestedCalculation`: the arguments of the call |
| `molecules` | The input `MoleculeCollection` |
| `rejections` | Tuple of `Rejection`: methods and parameter sets that were not used, with reasons |
| `input_warnings` | Tuple of `InputIssue`: cautions about the input molecules |
| `charge_warnings` | Tuple of `ChargeIssue`: cautions about the calculated charges |
| `failure_message` | Text of a numerical failure; otherwise `None` |
| `timings` | `CalculationTimings` |

(selecting-an-assignment)=
#### Selecting an assignment

`result.assignment(molecule=..., conformer=None)` returns one set of charges. `molecule` is required.
Omit `conformer` for charges that do not depend on coordinates; otherwise supply its position in the
input molecule, counted from 0.

| Exception | Cause |
| --- | --- |
| `IndexError` | The molecule index is out of range, the molecule has no charges, or a supplied conformer index is negative or absent from the per-conformer results |
| `ValueError` | You omitted `conformer` for charges calculated per conformer, or supplied it for charges that do not depend on coordinates |
| `TypeError` | An index is not an integer; booleans are not accepted |

### ChargeAssignment

| Attribute | Content |
| --- | --- |
| `values` | NumPy `float64` array of charges in elementary charges, in atom order |
| `molecule_index` | Index of the input molecule |
| `conformer_index` | Index of the conformer; `None` for methods that do not use coordinates |
| `atom_ids` | Atom IDs, as `Molecule.atom_ids` |
| `source` | `SourceIdentity` of the molecule |

`numpy.asarray(assignment)` returns `values`.

### ExecutedPlan

| Attribute | Content |
| --- | --- |
| `method` | `Method` |
| `parameter_set` | `ParameterSet`; `None` for methods without parameter sets |
| `options` | Values of all method options, including defaults |
| `execution` | `Execution` |
| `warnings` | Tuple of `ExecutionIssue` |
| `fixed_ions` | `FixedIons`; `None` without fixed ions |
| `settings` | `CalculationSettings` of the plan |

### Execution

| Attribute | Content |
| --- | --- |
| `mode` | `"full"`, `"cutoff"`, or `"cover"` |
| `radius` | Fragment radius in ångströms; `None` for full execution |

### CalculationSettings

The choices of a plan, which `settings` in `calculate()` and `assess()` apply to other molecules; see
[One model for all molecules](../interfaces/python.md#one-model-for-all-molecules). You get them
from `plan.settings` or `result.plan.settings`. Two settings are equal when all their choices are.

| Attribute | Content |
| --- | --- |
| `method` | `Method` |
| `parameter_set` | `ParameterSet`; `None` for methods without parameter sets |
| `options` | Values of all method options, including defaults |
| `parameter_matching` | `"strict"` or `"permissive"` |
| `execution` | `Execution` |
| `fixed_ions` | Tuple of the residue names requested for fixed ions, also those not found in the molecules |

### FixedIons and FixedAtomCharge

`FixedIons.sources` is a tuple of `FixedAtomCharge`, one per fixed ion, each with `component_id` (the
residue name it was selected by, such as `"MG"`), `molecule_index`, `atom_index`, and `charge`.

### CalculationTimings

`applicability_seconds` is the time spent checking methods and parameter sets, and
`computation_seconds` the time spent calculating charges.

### RequestedCalculation

The arguments of a call as ChargeFW read them: `method` as a `Method` and `parameter_set` as a
`ParameterSet`, or `None` when not given, `method_scope`, `options_by_method`, `parameter_matching`,
`execution`, `radius`, `cutoff_threshold`, `cover_threshold`, `threads`, and `fixed_ions` as a tuple
of residue names.

## Assessment

### Assessment

| Attribute | Content |
| --- | --- |
| `plans` | Tuple of `Plan` that can run, in selection order |
| `default_plan` | The `Plan` that `calculate()` would use; `None` when none can run |
| `rejections` | Tuple of `Rejection` |
| `input_warnings` | Tuple of `InputIssue` |
| `seconds` | Time taken |

### Plan

| Attribute | Content |
| --- | --- |
| `method`, `parameter_set`, `options`, `execution`, `warnings`, `settings` | As in `ExecutedPlan` |
| `estimate(*, threads=None)` | `ResourceEstimate` with `memory_bytes` and `operations`; see [Threads and estimates](../guide/large-structures.md#threads-and-estimates) |
| `calculate(*, threads=None, observer=None)` | `CalculationResult` for the molecules the plan was assessed with; raises and warns as `chargefw.calculate()` does |

### Rejection

| Attribute | Content |
| --- | --- |
| `method` | `Method` that was not used |
| `parameter_set` | `ParameterSet` that was not used; `None` when the method was rejected regardless of the parameter set |
| `execution` | `Execution` that was rejected; `None` when the method or parameter set was rejected before an execution mode was considered |
| `issues` | Tuple of `PrerequisiteIssue` and `ExecutionIssue`: the reasons |

### PrerequisiteIssue

A reason why a method or parameter set does not apply to the molecules.

| Attribute | Content |
| --- | --- |
| `kind` | One of the kinds below |
| `message` | Description, such as `molecule 1: parameter set 'EEM_Baek1991' has no atom parameter matching atom 1 (...)`; numbers in messages count from one |
| `molecule_index`, `atom_index`, `bond_index`, `conformer_index` | Where the problem is; `None` where it does not apply |
| `atom_indices`, `bond_indices` | Read-only NumPy `int64` arrays of every atom and bond the problem concerns; `atom_index` and `bond_index` are the ones the message names |

| Kind | Meaning |
| --- | --- |
| `parameter_classification_failed` | The parameter set has no entry matching some atoms or bonds of the molecule; one issue lists all such atoms, another all such bonds. [Permissive matching](parameter-format.md#strict-and-permissive-matching) may help |
| `missing_parameters` | The parameter set lacks values the method needs for every molecule, or no parameter set was available |
| `unknown_parameters` | The parameter set has values the method does not use, such as a set for another version of the method |
| `missing_feature` | The method needs coordinates and the molecule has none |
| `invalid_geometry` | Coordinates that are not finite, or two atoms at the same position |
| `unsupported_molecule` | The method does not handle this molecule, for example a charged component where it needs neutral ones, or elements it has no values for; one issue lists all such atoms |
| `invalid_options` | A method option has a value the method does not accept |
| `unsupported_fixed_ions` | The method does not support fixed ions |

### ExecutionIssue

A reason why an execution mode was not used, or a warning about the one that was. `kind` is
`"resource_threshold_exceeded"` or `"unsupported_execution_mode"`; it also has a `message` and a
`molecule_index` (or `None`).

### InputIssue

A caution about an input molecule, such as hydrogens that were probably not added; ChargeFW still
calculates.

| Attribute | Content |
| --- | --- |
| `kind` | `"hydrogens_missing"`: the molecule contains carbon but no hydrogen atoms; `"atoms_too_close"`: atoms of a conformer overlap. See [Preparing input](../guide/input.md) |
| `message` | Description |
| `molecule_index` | Index of the input molecule |
| `conformer_index` | Index of the conformer the caution concerns; `None` when it concerns the whole molecule |
| `atom_index` | Index of an atom the caution points to, for `"atoms_too_close"` the second atom of the closest pair; otherwise `None` |
| `atom_indices` | Read-only NumPy `int64` array of every atom the caution concerns, for `"atoms_too_close"` each atom of a close pair; empty for cautions about the whole molecule |

`"atoms_too_close"` appears only in the results of methods that use coordinates, not in
`Assessment.input_warnings`.

### ChargeIssue

A caution about calculated charges: ChargeFW returns them, but they are probably wrong; see
[Results](../guide/results.md#charges-that-are-probably-wrong).

| Attribute | Content |
| --- | --- |
| `kind` | `"charge_out_of_range"`: charges are more than 6 e away from the formal charges of their atoms |
| `message` | Description |
| `molecule_index` | Index of the input molecule |
| `conformer_index` | Index of the conformer of the charges; `None` for methods that do not use coordinates |
| `atom_index` | Index of the atom the message names, the one furthest from its formal charge |
| `atom_indices` | Read-only NumPy `int64` array of every affected atom, in atom order |

### ChargeFWWarning, InputWarning, and ChargeWarning

Python warnings for cautions that can make the charges wrong; see
[Warnings](../interfaces/python.md#warnings). `ChargeFWWarning` derives from `UserWarning`, and
`InputWarning` and `ChargeWarning` from `ChargeFWWarning`. `InputWarning` is raised for the
`ImportIssue` kinds `"unbonded_residue_atoms"` and `"symmetry_connections_ignored"` and for every
`InputIssue`; `ChargeWarning` for every `ChargeIssue`.

## Methods and parameter sets

### Method

| Attribute | Content |
| --- | --- |
| `id`, `name`, `full_name` | Such as `"peoe"`, `"PEOE"`, and `"Partial Equalization of Orbital Electronegativity"` |
| `publication` | DOI of the method's publication, or `None` |
| `notes` | Remarks on how ChargeFW implements the method |
| `recommended` | `True` for methods whose charges are recommended for general use; automatic selection uses only these unless `method_scope="all"` |
| `requires_coordinates` | `True` for methods that use 3D coordinates |
| `requires_bonds` | `True` for methods that use the bonds of the molecule |
| `time_complexity`, `memory_complexity` | Cost for $n$ atoms and $m$ bonds, such as `"O(n^3)"` |
| `supports_cutoff`, `supports_cover` | Whether the method can run the [approximations for large structures](../guide/large-structures.md) |
| `supports_fixed_ions` | Whether the method accepts [fixed ions](../guide/fixed-ions.md) |
| `options` | Mapping of option ID to `MethodOption` |
| `parameter_sets` | Mapping of ID to `ParameterSet` for this method, in the order automatic selection tries them |

### MethodOption

| Attribute | Content |
| --- | --- |
| `id`, `description` | Such as `"iters"` and its description |
| `type` | `"boolean"`, `"integer"`, `"floating_point"`, or `"string"` |
| `default` | Default value |
| `choices` | Allowed values for options with a fixed set of values; otherwise empty |
| `minimum`, `maximum` | Bounds, or `None` |
| `minimum_inclusive`, `maximum_inclusive` | Whether the bound itself is allowed |

### ParameterSet

| Attribute | Content |
| --- | --- |
| `id`, `name` | Such as `"EEM_Cheminf_b3lyp_npa"` |
| `method` | `Method` the set is for |
| `publication`, `notes` | Source of the values and remarks |
| `atom_keys` | Tuple of `AtomParameterKey`, each with `element`, `classifier`, and `type` |
| `bond_keys` | Tuple of `BondParameterKey`, each with `atoms` (two `AtomParameterKey`), `classifier`, and `type` |

The [parameter-set format](parameter-format.md#classifiers) explains classifiers and types.

### Element

`atomic_number`, `symbol`, `name`, `period`, and `group`, as described in
[Periodic table](../interfaces/python.md#periodic-table).

## Fitting

`chargefw.parameterization` fits parameter values to reference charges; [Fitting
parameters](../interfaces/fitting.md) shows how. It is experimental and may change in any release.

| Function | Returns |
| --- | --- |
| `make_parameter_template(molecules, method, *, atom_types, bond_types="plain", id="template")` | `ParameterSet` with an entry for every atom and bond type of the molecules, all values zero |
| `fit(evaluator, references, *, weights=None, start=None, free=None, bounds=None, profile="default", max_evaluations=None, observer=None)` | `FitResult` |

`references` and `weights` hold one array per assignment or one array for all charges; a reference of
`nan` marks a charge without one. `start` is a vector in layout order, a mapping from (key, parameter
name) to values, or a `ParameterSet` whose values are copied to the slots with the same key and
parameter name.

`chargefw.parameterization.profiles.PROFILES` maps a method ID to its `Profile` of default bounds,
start values, and fixed parameters: `"eem"`, `"sqe"`, `"sqeq0"`, and `"sqeqp"`. A profile's `bounds`,
`start`, and `terminal` are read-only. To customize a profile, use `dataclasses.replace()` and pass
the result as `fit(..., profile=...)`:

```python
from dataclasses import replace
from chargefw.parameterization.profiles import PROFILES

profile = replace(PROFILES["eem"], start={**PROFILES["eem"].start, "B": 3.0})
```

### Evaluator

`Evaluator(plan, *, threads=None)` for a `Plan` with full execution and a parameter set.

| Attribute | Content |
| --- | --- |
| `layout` | Tuple of `Slot`, the order of every value vector |
| `values` | The plan's parameter values |
| `matches` | Per slot, the number of atoms or bonds of its type; `None` for common parameters |
| `charges(values)` | The charges of every assignment, concatenated in assignment order |
| `jacobian(values, slots=None, bounds=None)` | Derivatives of the charges by finite differences, one column per slot |

A `Slot` has `scope` (`"common"`, `"atom"`, or `"bond"`), `key` (an `AtomParameterKey`, a
`BondParameterKey`, or `None`), and `name`.

### FitResult

`FitResult`, `Diagnostics`, and `FitProgress` compare and hash by object identity, not by their
contents. To compare fitted values, use `numpy.array_equal(first.values, second.values)` or
`numpy.allclose(first.values, second.values)`.

| Attribute | Content |
| --- | --- |
| `values`, `start` | Fitted and starting values of every slot |
| `free` | Indices of the fitted slots |
| `bounds` | Lower and upper values per slot |
| `success`, `message` | Whether and why the optimizer stopped |
| `evaluations`, `failed_evaluations` | Trial values the optimizer evaluated, and those whose calculation failed |
| `rmsd` | RMSD of the charges at `values` against the references, without weights, over the charges that have a reference |
| `diagnostics` | `Diagnostics` |
| `parameter_set(*, id, name="", publication="", notes="")` | The fitted values as a `ParameterSet` |

### Diagnostics and ChargeDiagnostics

| Attribute | Content |
| --- | --- |
| `r2` | R² of the charges against the references |
| `by_element` | Mapping of element symbol to `ChargeDiagnostics` |
| `by_atom_type` | Mapping of `AtomParameterKey` to `ChargeDiagnostics` |
| `matches` | Per free slot, the number of atoms or bonds of its type |
| `at_bounds` | Free slots at a bound |
| `singular_values`, `rank` | How strongly, and in how many independent combinations, the data determines the free slots |
| `weakest` | Free slots of the least determined combination |
| `failed_messages` | Messages of failed trial calculations, at most 20 |

`ChargeDiagnostics` has `count`, `rmsd`, `mean_error` (calculated minus reference), and `r2`. `r2`,
`by_element`, and `by_atom_type` cover only charges that have a reference, and `by_atom_type` only
types that such a charge belongs to; fixed ions belong to no type. `matches` counts every atom and
bond.

### FitObserver and FitProgress

Subclass `FitObserver` and override `on_evaluation(progress)` to receive a `FitProgress` after each
trial value of a fit, and `cancelled()` to return `True` when the fit should stop. See [Following a
long fit](../interfaces/fitting.md#following-a-long-fit).

| Attribute | Content |
| --- | --- |
| `evaluation`, `failed_evaluations` | Trial values evaluated so far, and those whose calculation failed |
| `objective`, `rmsd` | Weighted sum of squares and RMSD without weights at this trial value; `nan` when its calculation failed |
| `best_objective`, `best_rmsd`, `best_values` | The lowest weighted sum of squares so far, the RMSD there, and the values of every slot there |
| `elapsed_seconds` | Time since `fit()` was called |
| `parameter_set(*, id, name="", publication="", notes="")` | The best values so far as a `ParameterSet` |

## Progress

### CalculationObserver

Subclass it and override `on_progress(progress)` to receive progress, and `cancelled()` to return
`True` when the calculation should stop. See [Progress and
cancellation](../interfaces/python.md#progress-and-cancellation).

### CalculationProgress

| Attribute | Content |
| --- | --- |
| `phase` | `"computation_started"`, `"target_started"`, `"fragment_progress"`, `"target_finished"`, or `"computation_finished"` |
| `mode` | Execution mode |
| `method` | `Method` being calculated |
| `target_index`, `target_count` | Position of the current target and the number of targets; a target is a conformer for methods that use coordinates, otherwise a molecule |
| `completed_fragment_count`, `fragment_count` | Fragments done and in total, for cutoff and cover |
| `molecule_index`, `conformer_index` | The molecule and conformer being calculated |
| `elapsed_seconds` | Time since the calculation started |

## Exceptions

| Exception | Raised when |
| --- | --- |
| `chargefw.ChargeFWError` | Base class of the three below; `error.result` holds the `CalculationResult` |
| `chargefw.NoExecutablePlanError` | No method and parameter set apply with the requested settings |
| `chargefw.NumericalFailureError` | The calculation failed, for example did not converge |
| `chargefw.CalculationCancelledError` | The observer cancelled the calculation |
| `chargefw.parameterization.EvaluationError` | A calculation for parameter values failed during fitting or evaluation |

Invalid arguments raise `TypeError` or `ValueError`, files that cannot be read or written raise
`OSError`, and a calculation that runs out of memory raises `MemoryError`.
