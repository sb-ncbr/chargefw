# Python API reference

This page lists every public function, class, and attribute of the `chargefw` package.
[Python](../interfaces/python.md) shows how to use them together. All values are read-only, and all
indices of molecules, atoms, bonds, and conformers count from zero.

## Functions

| Function | Returns |
| --- | --- |
| `chargefw.calculate(molecules, *, method=None, parameter_set=None, ...)` | `CalculationResult`; the arguments are listed in [Calculating charges](../interfaces/python.md#calculating-charges) |
| `chargefw.calculate(molecules, plan, *, threads=None, observer=None)` | `CalculationResult` of a `Plan` from `assess()`; omit `plan` rather than passing `None` |
| `chargefw.assess(molecules, *, method=None, parameter_set=None, ...)` | `Assessment`; the same arguments as `calculate()` except `plan` and `observer` |
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
| `import_diagnostics` | Tuple of `ImportDiagnostic`: what the reader left out of the record |

### MoleculeCollection

`MoleculeCollection(molecules, name=None)` holds molecules in order. Index it and iterate over it
like a tuple; `molecules` returns the tuple. A [plan](#plan) works only with the collection it was
assessed with, which includes having the same `name`.

### ImportDiagnostic

| Attribute | Content |
| --- | --- |
| `code` | Kind of warning, such as `"partial_charges_ignored"` or `"coordinates_2d_ignored"` |
| `message` | Description |
| `line` | Line in the file, or `None` |

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
| `record_selection`, `conformer_selection`, `bond_strategy` | The `selection`, `conformers`, and `bonds` values used for reading, such as `"all"` and `"hybrid"`; `None` where the reader has no such option |
| `alternate_location_selection` | `"first-source-order"` for PDB and mmCIF, where the first [alternate location](formats.md#structural-selection) of each atom is read; otherwise `None` |

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
| `assignment(molecule=..., conformer=None)` | One assignment; see [Results](../interfaces/python.md#results) |
| `plan` | `ExecutedPlan` that was used, also after a numerical failure or cancellation; `None` when no plan could run |
| `requested` | `RequestedCalculation`: the arguments of the call |
| `molecules` | The input `MoleculeCollection` |
| `rejections` | Tuple of `Rejection`: methods and parameter sets that were not used, with reasons |
| `warnings` | Tuple of `ExecutionIssue`, such as full execution above the size limits |
| `failure_message` | Text of a numerical failure; otherwise `None` |
| `timings` | `CalculationTimings` |

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

### Execution

| Attribute | Content |
| --- | --- |
| `mode` | `"full"`, `"cutoff"`, or `"cover"` |
| `radius` | Fragment radius in ångströms; `None` for full execution |

### FixedIons and FixedAtomCharge

`FixedIons.sources` is a tuple of `FixedAtomCharge`, one per fixed ion, each with `component_id` (the
residue name it was selected by, such as `"MG"`), `molecule_index`, `atom_index`, and `charge`.

### CalculationTimings

`applicability_seconds` is the time spent checking methods and parameter sets, and
`computation_seconds` the time spent calculating charges.

### RequestedCalculation

The arguments of a call as ChargeFW read them: `method` and `parameter_set` as IDs or `None`,
`method_scope`, `options_by_method`, `parameter_matching`, `execution`, `radius`, `cutoff_threshold`,
`cover_threshold`, `threads`, and `fixed_ions` as a tuple of residue names.

## Assessment

### Assessment

| Attribute | Content |
| --- | --- |
| `plans` | Tuple of `Plan` that can run, in selection order |
| `default_plan` | The `Plan` that `calculate()` would use; `None` when none can run |
| `rejections` | Tuple of `Rejection` |
| `seconds` | Time taken |

### Plan

| Attribute | Content |
| --- | --- |
| `method`, `parameter_set`, `options`, `execution`, `warnings` | As in `ExecutedPlan` |
| `estimate(threads=None)` | `ResourceEstimate` with `memory_bytes` and `operations`; see [Threads and estimates](../guide/large-structures.md#threads-and-estimates) |

Pass a plan to `calculate()` after the molecules it was assessed with.

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

| Kind | Meaning |
| --- | --- |
| `parameter_classification_failed` | The parameter set has no entry matching an atom or bond; [permissive matching](parameter-format.md#strict-and-permissive-matching) may help |
| `missing_parameters` | The parameter set lacks values the method needs for every molecule, or no parameter set was available |
| `missing_feature` | The method needs coordinates and the molecule has none |
| `invalid_geometry` | Coordinates that are not finite, or two atoms at the same position |
| `unsupported_molecule` | The method does not handle this molecule, for example a charged component where it needs neutral ones |
| `invalid_options` | A method option has a value the method does not accept |
| `unsupported_fixed_ions` | The method does not support fixed ions |

### ExecutionIssue

A reason why an execution mode was not used, or a warning about the one that was. `kind` is
`"resource_threshold_exceeded"` or `"unsupported_execution_mode"`; it also has a `message` and a
`molecule_index` (or `None`).

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
| `method` | ID of the method the set is for |
| `publication`, `notes` | Source of the values and remarks |
| `atom_keys` | Tuple of `AtomParameterKey`, each with `element`, `classifier`, and `type` |
| `bond_keys` | Tuple of `BondParameterKey`, each with `atoms` (two `AtomParameterKey`), `classifier`, and `type` |

The [parameter-set format](parameter-format.md#classifiers) explains classifiers and types.

### Element

`atomic_number`, `symbol`, `name`, `period`, and `group`, as described in
[Periodic table](../interfaces/python.md#periodic-table).

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
| `method_id` | ID of the method |
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

Invalid arguments raise `TypeError` or `ValueError`, files that cannot be read or written raise
`OSError`, and a calculation that runs out of memory raises `MemoryError`.
