# Changelog

## 3.0.0rc5

This is the last release candidate for ChargeFW 3.0.0, the first stable release of the successor
of [ChargeFW2](https://github.com/sb-ncbr/ChargeFW2). ChargeFW has the same 20 methods and
parameter sets as ChargeFW2 in a new implementation: a Python package for CPython 3.12 and newer
on Linux, macOS, and Windows; a command line, also as the Docker image `sbncbr/chargefw`; and a C++
library. [Differences from ChargeFW2](https://chargefw.readthedocs.io/en/latest/about/chargefw2.html)
describes what is new and how to move over. The version number continues from ChargeFW2, so the
`chargefw` Python package moves from 0.1.7 to 3.0.0. Install the release candidate with
`pip install --pre chargefw`.

The changes since 3.0.0rc4 follow.

### Installation

- The wheels are tested on CPython 3.15. Gemmi 0.7.4 has no wheel for it, and building Gemmi from
  source fails with nanobind 3, so the `gemmi` extra needs a build constraint there;
  [Installation](https://chargefw.readthedocs.io/en/latest/getting-started/installation.html#python-package)
  shows the command.
- Native installations put the documentation, license notices, and JSON Schemas directly in
  `share/doc/chargefw/` instead of `share/doc/chargefw/chargefw/`.

### Command line

- `chargefw --version` prints the version, such as `ChargeFW 3.0.0`.

### Output

- Breaking change: the JSON Schemas are `result-1.schema.json` and `molecule-input-1.schema.json`,
  without the minor version, because each one accepts every 1.x document. Their `$id` is the address
  where the documentation publishes them, such as
  `https://chargefw.readthedocs.io/en/stable/schemas/result-1.schema.json`, instead of a
  `urn:chargefw:schema:` name.

## 3.0.0 pre-releases

Each pre-release lists its changes since the previous one. A breaking change in them breaks only
code written for the previous pre-release.

### 3.0.0rc4

#### Input

- Breaking change: the options for reading structures have the same names everywhere, those of
  `chargefw.io.read()`: `selection`, `bonds`, and `conformers`. On the command line,
  `--structural-selection` and `--structural-bonds` are now `--selection` and `--bonds`. In result
  JSON, `import.policy` has `selection`, `bonds`, `conformers`, and `alternate_locations` instead of
  `record_selection`, `bond_strategy`, `conformer_selection`, and `alternate_location_selection`. In
  Python, these four are in `SourceMapping.policy`, a new `chargefw.SourcePolicy`; in C++, in
  `MoleculeImportMetadata::policy`, an `io::ImportPolicy`, and the reader options'
  `bond_strategy` is `bonds` (`io::StructureOptions`, formerly `adapters::gemmi::InputOptions`).
- Breaking change: parameter-set files are read strictly. An unknown section, key member, or entry
  member is an error, and so is a parameter name the method does not use, reported as the new
  prerequisite issue kind `unknown_parameters`. A set made for another version could otherwise give
  charges without values it contains. `metadata` may still hold other members, and its new optional
  `format_version` must be `1`.
- Errors reading MOL, SDF, and MOL2 input name the line and record, and the source file when supplied.
- MOL, SDF, MOL2, PDB, and mmCIF input must be ASCII text. A character outside ASCII on any line,
  including comments and SDF data items, is an error that names the line.
- MOL and SDF records are read for their atoms, bonds, and charges only; everything else is ignored.
  V3000 `SGROUP`, `3D`, and other blocks and V2000 lines such as atom aliases, which were errors, are
  now ignored. Breaking change: the `v2000_property_ignored` warning was removed.
- Breaking change: molecule JSON input rejects members its format does not define, at every level,
  so a misspelled member such as `bond` no longer gives a molecule without bonds. The molecule input
  schema sets `additionalProperties: false` to match.

#### Installation

- A patch release, such as 3.0.1, contains fixes only; a minor release, such as 3.1.0, may also
  change documented interfaces, with every breaking change marked here. Pin `chargefw~=3.0.0` to stay
  within one minor version. `chargefw.parameterization` is experimental and may change in any release.
  [Versions and changes](https://chargefw.readthedocs.io/en/latest/about/versions.html) has the details.
- On x86-64, the wheels, the Docker image, and source builds target x86-64-v3 (AVX2 and FMA: Intel
  Haswell, AMD Excavator, or newer). Full calculations on large structures run about twice as fast
  with the wheels. On an older CPU, `import chargefw` raises an `ImportError`; build from source with
  `CHARGEFW_ENABLE_NATIVE_OPTIMIZATIONS=ON`, which now defaults to `OFF`.

#### Python

- Breaking change: `FixedAtomCharge` and `FixedIons` take their attributes as keyword arguments
  only, as the other value classes, for example
  `FixedAtomCharge(molecule_index=0, atom_index=3, charge=2.0, component_id="MG")`.
- Breaking change: `Method.publication` is `""` for a method without a publication, not `None`,
  as for parameter sets. `PeriodicTable` can no longer be constructed; use
  `chargefw.periodic_table`.
- Breaking change: `Molecule.import_issues` is now `import_warnings`, as `input_warnings`,
  `charge_warnings`, and `plan.warnings` name the other cautions.
- Breaking change: both `attach_charges()` functions, in `chargefw.io.gemmi` and
  `chargefw.io.rdkit`, take `if_exists="error"` (default) or `"replace"` instead of
  `overwrite=False` or `True`. The values are `chargefw.io.IfExists`.
- `chargefw.methods`, `chargefw.parameter_sets`, `Method.options`, and `Method.parameter_sets` are
  typed as `Mapping[str, …]`, as documented, instead of private catalog classes.
- Breaking change: `cutoff_threshold` and `cover_threshold` take `"unlimited"` for no limit, as
  on the command line and in result JSON, instead of `None`, which is now rejected.
  `RequestedCalculation` holds `"unlimited"` too.
- Breaking change: `SourceComponentInstance.atom_indices` is a read-only NumPy `int64` array, as
  the `atom_indices` of the issues, not a tuple, so it indexes charges directly. Components
  compare by identity.
- Breaking change: slicing a `MoleculeCollection` returns a `MoleculeCollection` with the same
  `name`, not a tuple. Its `name` argument is keyword-only, defaults to `""`, and no longer
  accepts `None`.
- Breaking change: `CalculationResult.assignment()` takes `molecule_index` and `conformer_index`
  instead of `molecule` and `conformer`, and `chargefw.io.rdkit.attach_charges()` takes
  `conformer_index` instead of `conformer`, as the attributes of `ChargeAssignment` are named.
- Breaking change: `Molecule()` takes only `atomic_numbers` by position; the other arguments are
  keyword-only. The record it came from is given only by `source_name`, `record_index`, and
  `record_id`; the `source` argument is removed. `name` and `source_name` default to `""` and
  no longer accept `None`. `SourceIdentity` has only `source_name`, as the argument of
  `Molecule()` and `chargefw.io.parse()`; its `source` field is removed.
- Breaking change: `chargefw.ChargeFWError` is the base class of every ChargeFW error and no longer
  a `RuntimeError`. A failed calculation raises a `chargefw.CalculationError`, which holds `result`
  and is the base class of `NoExecutablePlanError`, `NumericalFailureError`, and
  `CalculationCancelledError`. Molecules or a parameter set that `io.read()`, `io.parse()`, or
  `read_parameter_set()` cannot read raise the new `chargefw.ParseError`, also a `ValueError`,
  instead of `RuntimeError`, `ValueError`, or `IndexError`. `read_parameter_set()` raises
  `OSError` for a file it cannot open. `parameterization.EvaluationError` derives from
  `ChargeFWError` instead of `RuntimeError`.
- Breaking change: the modules `chargefw.core`, `chargefw.charges`, and `chargefw.calculation` are
  private; import their classes and functions from `chargefw`, as the documentation does. Every class
  and function of `chargefw` now reports `chargefw` as its module. `chargefw.io` and
  `chargefw.parameterization` stay public.
- Breaking change: RDKit `attach_charges()` delegates explicit conformer selection to
  `CalculationResult.assignment()`. Negative or missing conformers now raise `IndexError` instead of
  `ValueError`; explicit selection no longer falls back to geometry-independent charges and raises
  `ValueError` instead. Omitting `conformer` still selects a single assignment automatically,
  including a conformer-specific assignment.
- Breaking change: the value classes take their attributes as keyword arguments only, for example
  `SourceIdentity(source_name="input.sdf", record_index=0)` and
  `AtomParameterKey(element="C", classifier="bonded", type="COO")`. This covers the issue classes,
  `SourceIdentity` and the other `Source*` classes, `ChargeAssignment`, `ExecutedPlan`, `Execution`,
  `ResourceEstimate`, `CalculationProgress`, `CalculationTimings`, `Rejection`, `Method`,
  `MethodOption`, `ParameterSet`, `AtomParameterKey`, `BondParameterKey`, and `Element`. Their
  attributes no longer match positional patterns in `match` statements. `Plan.estimate()` takes
  `threads` as a keyword argument only.
- Breaking change: an attribute named `method` is always a `Method`, and one named `parameter_set`
  always a `ParameterSet`. `ParameterSet.method`, `RequestedCalculation.method`, and
  `RequestedCalculation.parameter_set` were IDs, and `CalculationProgress.method_id` is now
  `CalculationProgress.method`; use `.id` for the ID. `read_parameter_set()` rejects a file for a
  method ChargeFW does not have; such a set could be read before but never calculated.
- Breaking change: `ParameterSet`, `CalculationSettings`, `Plan`, `Assessment`, and
  `CalculationResult` cannot be constructed; calling the class raises `TypeError`. Get them from
  `chargefw.parameter_sets`, `read_parameter_set()`, `assess()`, and `calculate()`. `ParameterSet` is
  no longer a dataclass, so `dataclasses.replace()` no longer makes a copy that names a parameter set
  that does not exist. Two `CalculationSettings` with the same choices are now equal.
- Breaking change: `chargefw.calculate()` no longer takes a plan. Calculate a plan from `assess()` with
  `plan.calculate(threads=None, observer=None)`; it calculates the molecules the plan was assessed with.
- `chargefw.parameterization` starts support for fitting parameter values: `make_parameter_template()`
  makes a parameter set with an entry for every atom and bond type of your molecules, `Evaluator`
  calculates the charges of an assessed plan, and their derivatives, for any vector of values, and
  `fit()` fits the values to reference charges with SciPy's least squares, using exact derivatives
  for EEM, SQE, SQE+q0, and SQE+qp, whose cost then does not grow with the number of charges in each
  step of the optimizer; the result's
  `parameter_set()` turns the fitted values into a parameter set. Profiles for EEM and the SQE
  methods give default bounds, starting values, and fixed parameters; the SQE profiles fit each
  hardness as √π·hardness·width − 1 ≥ 0, which keeps the system stable for every molecule and
  geometry, and start the electronegativities (EEM's A, and SQE+qp's q0) from their best values for
  the start of the other parameters, so they need no start. Fit diagnostics report R² and the error
  per element, parameters at a bound, failed calculations, which parameters the data cannot
  determine, and how many atoms or bonds each type has. Install SciPy 1.18 or newer with the `parameterization` extra.
- `chargefw.read_parameter_set()` reads a parameter set from your own JSON file, which `assess()` and
  `calculate()` accept as `parameter_set`. Its ID must differ from the IDs of the bundled sets.
  `chargefw.write_parameter_set()` writes any parameter set to a JSON file.

#### C++

- Breaking change: the readers, writers, records, and results are in `chargefw::io`, as Python's
  `chargefw.io`, instead of `chargefw::adapters` and its per-format namespaces: `io::SdfReader`
  instead of `adapters::native::sdf_input::SdfReader`, `io::MmcifWriter` instead of
  `adapters::gemmi::mmcif_output::MmcifWriter`, and so on. Their headers are in `chargefw/io/`
  without the `native/` and `gemmi/` subdirectories, for example `<chargefw/io/sdf_input.h>`, and one
  `<chargefw/io/all.h>` replaces the two umbrella headers. `adapters::gemmi::InputOptions` is
  `io::StructureOptions`, in `<chargefw/io/structure_options.h>` with `RecordSelection`,
  `BondStrategy`, and their conversion functions. Public names no longer hide the `::gemmi`
  namespace.
- `io::make_requested_provenance()` creates the `RequestedCalculationProvenance` of a request, with
  the thread limit and the residue names of the fixed ions, instead of copying its fields by hand.
- Breaking change: `ImportedMoleculeRecord::issues` is now `import_warnings`, as in Python.
- Breaking change: `ExecutionMetrics::peak_resident_memory_mb` is now `peak_resident_memory_bytes`,
  a `std::uint64_t`.
- Breaking change: `MoleculeRecordIdentity::source` is now `source_name`, as in Python and result
  JSON.
- The MOL, SDF, and MOL2 readers throw `std::runtime_error` for an unknown element symbol, as
  documented for input they cannot parse, instead of `std::out_of_range`.
- Breaking change: the request and its provenance use the names of the Python arguments.
  `AssessmentRequest::classification_options` is `parameter_matching`, a new
  `parameters::ParameterMatching` (`strict` or `permissive`), which also replaces
  `ClassificationOptions::permissive_types`. `AssessmentRequest::execution_selection` is
  `execution`, `ExecutionSelection::kind()` is `mode()`, and `ExecutionSelectionKind` is
  `ExecutionSelectionMode`. `ResourcePolicy` has `cutoff_threshold` and `cover_threshold`, and the
  defaults are `default_cutoff_threshold` and `default_cover_threshold`. In
  `RequestedCalculationProvenance`, `permissive_types`, `max_threads`, `execution_kind`, and
  `execution_radius` are `parameter_matching`, `threads`, `execution`, and `radius`.
  `ResultStatus::invalid_request` is `invalid_input_or_request`, as result JSON writes it.
- `find_package(chargefw 3.0)` accepts only 3.0.x, because a minor release may change interfaces.
- Breaking change: the CMake target is `chargefw::chargefw` instead of `chargefw::core`, and the
  library file is `libchargefw.so` instead of `libchargefw_core.so`; `chargefw::core` named both the
  whole library and the namespace of atoms and molecules, which keeps its name. Link
  `chargefw::chargefw` after `find_package(chargefw)`.
- The installed C++ library no longer requires nlohmann/json or Gemmi development packages in downstream
  projects. JSON and mmCIF reader state is private to the library; parsing behavior is unchanged.
  Breaking change: the library target no longer supplies these dependencies' headers or CMake targets.
  Applications that use either dependency directly must find and link it explicitly. Rebuild C++
  consumers because the reader layouts changed. Fetched native installations retain the Gemmi and
  oneTBB runtime libraries but no longer install third-party development packages.
- The shared library carries its full release version in its name, for example
  `libchargefw.so.3.0.0rc3`. A program needs the library of the release it was built against
  and fails to start without it, instead of loading another release. Programs linked against earlier
  releases ask for the unversioned name and get no such protection; rebuild them, and rebuild C++
  programs against every release.
- Breaking change: `<chargefw/config.h>` is replaced by `<chargefw/version.h>`, which defines
  `CHARGEFW_VERSION_STRING` as before and the numbers `CHARGEFW_VERSION_MAJOR`,
  `CHARGEFW_VERSION_MINOR`, and `CHARGEFW_VERSION_PATCH`. The parameter-directory macros are no
  longer installed.
- Breaking change: removed `<chargefw/methods/method_applicability.h>` with
  `methods::ApplicableMethod`, `methods::RejectedCandidate`, and `methods::ApplicabilityResult`. No
  public function produced them. Read the plans and rejections of an assessment instead.
- Breaking change: removed `features::TopologyFeatures::adjacency()`. Use
  `neighbor_indices(atom_index)` for each atom.
- Breaking change: removed `methods::validate_method_option_schema()` and
  `parameters::validate_parameter_classification()`. `make_default_options()` and
  `validate_method_options()` still reject an invalid option schema, and `CalculationInput` still
  rejects a classification that does not fit its molecule and parameter set.
  `<chargefw/parameters/classification/parameter_classification.h>` no longer includes the molecule
  and parameter-set headers.
- Breaking change: three names now match what they name.
  `calculation::EffectiveCalculation::execution_issues` is `warnings`, as in
  `ExecutionPlan::warnings()`; `<chargefw/charges/charge_collection.h>` is
  `<chargefw/charges/charge_set.h>`, after the `ChargeSet` it declares; and
  `methods::complexity_notation(term)` is `methods::to_string(term)`.
- Enumerations and `std::variant` types, such as `calculation::RejectionIssue`, may gain values and
  alternatives in later releases; give a `switch` a `default` and `std::visit` a fallback overload.
- `methods::Method` is the interface of the built-in methods, and only ChargeFW implements it:
  requests select methods by ID from `methods::method_registry()`, so a class derived elsewhere is
  never selected. Later releases may add pure virtual functions to it.
- `methods::MethodRegistry::require(id)` returns a `const Method&` or throws
  `std::invalid_argument` if the method is not registered. `find(id)` remains a nullable lookup.
- Breaking change: `AtomParameterAccessor::operator[]` and `BondParameterAccessor::operator[]`
  now require a valid atom or bond index, matching the other indexed containers. Use `at()` for
  bounds checking and `std::out_of_range` on an invalid index. Named parameter lookup remains
  checked in both forms.
- `methods::CalculationInput` now rejects supplied parameter classifications with mismatched
  atom or bond counts or references to nonexistent parameter entries.
- Breaking change: removed `calculation::CalculationCancelled` from the public headers. Public
  `calculate()` overloads report cooperative cancellation as `ExecutionStatus::cancelled`; check
  the returned status instead of catching this internal exception.
- Breaking change: `methods::ExecutionIssue` and `methods::ExecutionIssueKind` are now
  `calculation::ExecutionIssue` and `calculation::ExecutionIssueKind`, declared in
  `<chargefw/calculation/execution_issues.h>`. Removed `methods::ExecutionAvailability`,
  `methods::ExecutionAssessment`, and `methods::ApplicableMethod::execution_assessments`; read the
  plans and rejections of an assessment instead. Rebuild downstream C++ code.
- Breaking change: an `ExecutionPlan` keeps the molecules and parameter sets of its assessment, so
  the functions that took an assessment and one of its plans take only the plan:
  `calculation::calculate(plan, ...)`, `calculate(plan, candidate, ...)`, `estimate(plan, ...)`,
  `parameterization::normal_equations(plan, ...)`, `entry_counts(plan)`, `atom_entries(plan)`, and
  `single_bond_entries(plan)`. Drop the assessment argument. A copy of a plan stays valid after its
  assessment is destroyed. `plans()` and `default_plan()` no longer compile on a temporary
  `AssessmentResult`; store the assessment first. Rebuild downstream C++ code.
- Breaking change: `features::PreparedMolecule` and `features::TopologyFeatures` can be neither
  copied nor moved, and `features::PreparedMoleculeCollection` can be moved but not copied. Build
  them where they are used, and pass references. Rebuild downstream C++ code.
- Breaking change: removed `PdbReader::options()` and `MmcifReader::options()`. Keep the
  `StructureOptions` supplied to the constructor if you need to inspect them later.
- Breaking change: `parameters::ClassificationResult` is no longer default-constructible. Obtain a
  result from `try_classify_parameters()`, or construct it explicitly from a classification or issues.
- Breaking change: `methods::PrerequisiteResult` now exposes a public `issues` vector. Replace
  `issues()` with `issues` and `add(issue)` with `issues.push_back(issue)`. Its boolean conversion
  still reports whether the vector is empty. Rebuild downstream C++ code.
- Breaking change: removed `parameters::ParameterSet::publication()` and `notes()`; use
  `metadata().publication` and `metadata().notes`. Removed `calculation::ExecutionResult::cancelled()`;
  compare `status` with `calculation::ExecutionStatus::cancelled` instead. Rebuild downstream C++ code.
- Breaking change: removed `features::TopologyFeatures::bond_index_between()` and `are_bonded()`.
  Use `neighbor_indices()` to check whether two atoms are bonded. To find the bond between them,
  look among `incident_bond_indices()` of one atom for the bond whose other atom is the second.
  Rebuild downstream C++ code.
- Breaking change: removed `methods::CalculationInput::has_geometry()`, `geometry_if_available()`,
  `has_parameters()`, and `parameters_if_available()`. Use `geometry()` and `parameters()` to access
  supplied values; these return references and throw `std::logic_error` when the value is missing.
  Rebuild downstream C++ code.
- Breaking source and ABI change: `RequestedCalculationProvenance::method_scope` and `execution`
  (formerly `execution_kind`) use `calculation::MethodScope` and
  `calculation::ExecutionSelectionMode` instead of strings. Replace string assignments and
  comparisons with enum values, using `ExecutionSelectionMode::automatic` for `"auto"`, and rebuild
  downstream C++ code. The defaults
  remain recommended methods and automatic execution; result JSON and Python strings are unchanged.
- Breaking change: the enumerators of `parameters::ClassificationIssueKind`,
  `parameters::AtomParameterClassificationKind`, and `parameters::BondParameterClassificationKind`
  are lower case, as in the other public enums: `PLAIN` is now `plain`, `HIGHEST_BOND_ORDER`
  `highest_bond_order`, and so on. Rename them in calls and rebuild downstream C++ code.
- Breaking change: removed `AtomParameters::entry()`, `BondParameters::entry()`,
  `AtomParameterClassification::parameter_entry_index()`, and
  `BondParameterClassification::parameter_entry_index()` from `chargefw::parameters`.
  Replace calls with `at()`, which retains the same bounds checking, and rebuild downstream C++ code.
- Breaking change: the single-molecule `methods::check_parameter_prerequisites()` overload and its
  `ParameterPrerequisiteInput` and `ParameterPrerequisiteResult` types were removed from the public
  headers. Use the collection overload with a `PreparedMoleculeCollection` and rebuild downstream
  C++ code.
- Breaking change: `parameters::classify_parameters()` and `parameters::try_classify_parameters()`
  derive the molecule from `TopologyFeatures`. Remove the first (`molecule`) argument from calls
  and rebuild downstream C++ code.
- Breaking change: `JsonWriter::write`, `Mol2Writer::write`, `PqrWriter::write`, and both
  `MmcifWriter::write` overloads no longer accept `generator_version`. Writers always record ChargeFW
  and the compiled library version. Pass optional JSON execution metrics as the second argument.
- The new `chargefw::parameterization` namespace has what fitting parameter values needs:
  `make_parameter_template()` makes a parameter set with an entry for every atom and bond type of a
  molecule collection under chosen classifiers, with the given parameter names and zero values;
  `parameter_layout()` and `parameter_values()` list the values of a set in a fixed order, and
  `with_parameter_values()` copies a set with new values; `normal_equations(assessment, plan,
  candidate, references, weights)` returns the charges, the weighted sum of squares, Jᵀ diag(w²) J,
  and Jᵀ diag(w²) r of a least-squares fit of the charges, from exact derivatives, for the methods
  `has_exact_derivatives()` names: EEM, SQE, SQE+q0, and SQE+qp; `entry_counts()` counts the atoms
  and bonds of a plan that each entry of its set matches, and `single_bond_entries()` reports the
  atom entries whose atoms all have one bond, with the bond entries of those bonds.
- `parameters::compatible()` tells whether two parameter sets differ only in their values and
  descriptive metadata, and parameter keys can be compared with `==`.
- `calculation::to_string()` names a `CalculationPhase`.
- `io::to_string()` names a `MolecularSourceFormat` and a `SourceConnectivity`, and
  `molecular_source_format_from_string()` reads a format name.
- `calculation::calculate(assessment, plan, candidate)` runs a plan with other values for its parameter
  set, without checking the molecules again.
- `parameters::save_parameter_set_json()` writes a parameter set in the JSON format that the loading
  functions read.
- `MethodRegistry` cannot be constructed; `method_registry()` is the catalog of built-in methods.
  `find_applicable_methods()`, `ApplicabilityRequest`, `CalculationRequest`, `CalculationResult`, and
  the `calculate()` overload that takes a `CalculationRequest` are removed. Calculate charges with
  `assess()` and `calculate()`.

#### Calculation

- Breaking change: method options belong to the method you name, and are called method options
  everywhere. In Python, `method_options={"iters": 8}` with `method=` replaces `options` and
  `options_by_method`, and `Plan`, `ExecutedPlan`, and `CalculationSettings` have `method_options`
  instead of `options`. On the command line, `--method-option iters=8` needs `--method` instead of
  naming the method as in `peoe.iters=8`. In C++, `AssessmentRequest::method_options` and
  `RequestedCalculationProvenance::method_options` are one `MethodOptions`. In result JSON,
  `requested.method_options` and `effective.method_options` map option IDs to values instead of
  method IDs to options. Options for automatic selection are no longer accepted: they applied only
  if that method happened to be chosen.
- Breaking change: automatic selection of an EEM parameter set follows ChargeFW2 again. It tries the
  NEEMP sets first, then Cheminf, Bult, Ouy, Svob, and Baek, and the Ionescu sets last; earlier 3.0
  pre-releases tried Ionescu first. `--method eem` without a parameter set therefore now uses
  `EEM_NEEMP_ccd2016_npa` for molecules it covers. Name the parameter set, such as
  `EEM_Ionescu2013_mpa_gas`, to keep the charges of earlier pre-releases.
- SQE, SQE+q0, and SQE+qp calculations fail with `numerical_failure` when the split-charge system
  cannot be factorized as positive definite, as when the energy has no unique minimum at the
  molecule's geometry. Earlier versions returned the solution without warning, with charges that could
  be far too large. In cutoff and cover execution, each fragment is checked. Solving the system takes
  about half the time.

#### Command line

- Breaking change: `--fixed-ions` takes one residue name per occurrence, as `--method-option`
  takes one option: write `--fixed-ions MG --fixed-ions CA`, not `--fixed-ions MG CA`.
- Breaking change: options are named as the arguments of `chargefw.calculate()` in Python.
  `--permissive-types` is now `--parameter-matching permissive` (default `strict`), and
  `--cutoff-atom-threshold` and `--cover-atom-threshold` are now `--cutoff-threshold` and
  `--cover-threshold`.
- Breaking change: a result JSON or output directory that cannot be written exits with status 6, as
  MOL2, mmCIF, and PQR do, instead of 2. Status 2 now means only invalid input or options.

#### Output

- Breaking change: result JSON reports `execution_metrics.peak_resident_memory_bytes`, an integer,
  instead of `peak_resident_memory_mb`, which was in MiB.
- Breaking change: a record's `input.source` in result JSON is now `input.source_name`, as in
  Python.
- Breaking change: result JSON has no `calculation_provenance.effective.warnings`. It repeated
  the messages of the plan's warnings, which are diagnostics with their codes.
- Breaking change: `calculation_provenance.requested` in result JSON is flat and uses the names of
  the arguments of `chargefw.calculate()`: `parameter_matching` (`"strict"` or `"permissive"`)
  replaces `classification.permissive_types`; `execution` is a string and `radius` a number or
  `null` beside it, instead of `execution.kind` and `execution.radius_angstrom`; `cutoff_threshold`,
  `cover_threshold`, and `threads` replace the `resource_policy` object; and the new `fixed_ions`
  lists the residue names requested as fixed ions. `effective.execution.radius_angstrom` is
  `radius`.
- The result JSON Schema accepts every 1.x document: `schema_version` is any `1.x`, objects may have
  members it does not list, and status, severity, input format, selections, execution kind and mode,
  and method scope may have values it does not list; assignment scope and source connectivity stay
  fixed. Later 1.x versions add members and values but never remove or rename them;
  [File formats](https://chargefw.readthedocs.io/en/latest/reference/formats.html#chargefw-result-json-10)
  says how to read them. Documents written by this version are unchanged.

### 3.0.0rc3

#### Input

- Histidine and arginine get bond orders that fit their hydrogens and formal charges, using the
  structures of the wwPDB protonation variants. Before, a neutral histidine protonated on ND1 got the
  bonds of the cation and wrong charges.
- Nucleotide templates connect the sugar atoms. Before, PDB nucleic acids and mmCIF files without
  `_chem_comp_bond` had their bases and phosphates as separate fragments, and automatic selection fell
  back to methods that do not use bonds.
- The N-terminal hydrogens `H1` and `H3` are bonded only in amino acids without a peptide link to the
  previous residue. Before, a hydrogen that Hydride names `H1` on the carboxyl carbon of a C-terminal
  residue was bonded to the backbone nitrogen.
- MOL2 files are read with the formal charges in their `UNITY_ATOM_ATTR` section, which Open Babel
  writes for charged atoms, and with their sections in any order. Before, such files were refused.
- Records with aromatic bonds, which are read as single bonds, get an `aromatic_bonds_converted`
  warning.
- Reading takes less memory, especially for many small molecules and for PDB files with many models
  read with `--conformers first`. Bonds stated in large PDB and mmCIF files are read in linear time.

#### Calculation

- A molecule with carbon but no hydrogen atoms gets a `hydrogens_missing` warning, because its
  hydrogens were most likely not added.
- With a method that uses coordinates, a conformer whose atoms are closer than 0.7 times the sum of
  their covalent radii gets an `atoms_too_close` warning, which names the pair closest relative to
  that sum and lists every atom of a close pair.
- A calculated charge more than 6 e away from its atom's formal charge gives a
  `charge_out_of_range` warning, one per molecule or conformer, which names the furthest atom. Fixed
  ions are not checked.
- A parameter set without values for some atoms or bonds of a molecule is rejected with one issue for
  the atoms, which counts them by element, and one for the bonds, each listing all of them, instead of
  one issue per atom or bond. A method that does not support some elements gives one issue for all such
  atoms. A crystal structure without hydrogens gave tens of thousands of issues; it now gives a few.
- Fixed ions in large systems are checked against the other atoms much faster.
- QEq uses the Nishimoto-Mataga Coulomb term by default instead of DasGupta-Huzinaga, so QEq charges
  change. Set `overlap_term` to `DasGupta-Huzinaga` for the previous results.
- The QEq `overlap_term` choice `Louwen-Vogt` is removed. It used the hardness in place of the
  Coulomb integral at zero distance, which made charges diverge.
- SQE+qp normalizes its parameterized charges in each connected component, so that every protein
  chain, ligand, and water molecule keeps its formal charge. Before, one normalization over the whole
  structure gave each component a total that depended on the others: in a crystal structure, each
  water got about +0.07 e and protein chains were off by several electrons, and a ligand with charge −2
  came out near −0.1. Results for molecules with one component do not change.

#### Output

- MOL2 output writes SYBYL atom types, such as `C.3`, `N.am`, or `O.co2`, derived from the bonds and
  formal charges, and the formal charges in a `UNITY_ATOM_ATTR` section that Open Babel, RDKit, and
  ChargeFW read. Before, atom types were element symbols and formal charges were lost. Aromatic rings
  are written in their Kekulé form.
- mmCIF output writes structures whose histidines or arginines are in different protonation states.
  The bonds whose order differs are written with order `delo`, and the bond rows of a residue name
  follow its template. Reading the file back gives the same bonds.
- The warnings of the executed plan also appear in result JSON as warning diagnostics with their
  code, such as `resource_threshold_exceeded`: in the record of the molecule they concern, or in the
  document's diagnostics. `effective.warnings` still lists them as text.
- Writing results of rejected large structures takes linear instead of quadratic time.
- mmCIF output accepts charges from −10 to +10, the range of version 1.1 of the charge dictionary.
  Before, charges outside [−5, 5] were refused, which excluded fixed ions such as W⁶⁺.
- A result JSON diagnostic about several atoms or bonds lists all of them in `atom_indices` or
  `bond_indices`.
- For PDB and mmCIF input, result JSON and command-line messages name the residue and chain of the
  atoms a diagnostic points to, such as `atoms 3707 and 3708 are in PRO D 58`.

#### Python

- Problems that can make the charges wrong raise Python warnings: `chargefw.InputWarning` for
  missing hydrogens, atoms too close, and missing bonds in a structure file, and
  `chargefw.ChargeWarning` for charges far from the formal charges. One warning per kind and call
  counts the molecules. Before, they were only recorded in the result, where a script that did not
  look missed them.
- `molecule.import_diagnostics` is renamed to `molecule.import_issues`, and its values from
  `chargefw.ImportDiagnostic` to `chargefw.ImportIssue`, with a `kind` instead of a `code`. The kinds
  are the same strings, such as `"partial_charges_ignored"`. An `ImportIssue` lists the atoms and
  bonds it concerns in `atom_indices` and `bond_indices`.
- `CalculationResult.warnings` is removed; use `result.plan.warnings`. `result.plan` is `None` when no
  plan was recorded.
- `Assessment.input_warnings` and `CalculationResult.input_warnings` list `chargefw.InputIssue`
  values, such as missing hydrogens. An `InputIssue` about atoms has the `atom_index` its message
  names and the `atom_indices` it concerns.
- `CalculationResult.charge_warnings` lists `chargefw.ChargeIssue` values, cautions about the
  calculated charges, each with the `atom_indices` it concerns.

#### Command line

- An unexpected error during assessment is reported as an internal error, exit status 1. Before, it
  was reported as a numerical failure, exit status 4, with result JSON.
- Each warning of the executed plan, such as a full calculation above the size limit, is printed once,
  before the calculation. Before, it was printed once per molecule.
- A warning about a molecule names it, such as `Warning: molecule 1 ('1CRN'): …`. Several molecules
  with the same warning give one line with their number and the first of them. Before, the line did
  not say which molecule it concerned.

#### C++

- `adapters::MoleculeRecordDiagnostic` is replaced by `adapters::ImportIssue`, with an
  `ImportIssueKind` instead of a string code, and `ImportedMoleculeRecord::diagnostics` is renamed to
  `issues`.
- `ChargeCalculationResult::execution()` returns a pointer, which is null for a request that failed
  with an exception; create such a result with `make_failed_request_result()`. `status()` returns the
  new `ResultStatus`, and `charges()` the charges of a successful calculation.
- `calculation::ExecutionStatus::invalid_input_or_request` is removed; the library never returned it.
- `SourceAtomReference::structural_labels` is a `std::shared_ptr<const SourceStructuralLabels>`
  instead of a `std::optional`. Readers share one object among sites with equal labels; replace it
  instead of editing it.
- `calculation::InputIssue`, `AssessmentResult::input_warnings()`, and
  `ExecutionResult::input_warnings` report cautions about input molecules, such as missing hydrogens.
- `calculation::ChargeIssue` and `ExecutionResult::charge_warnings` report cautions about calculated
  charges. `features::CloseAtoms::atoms` lists the atoms of close pairs.
- `methods::PrerequisiteIssue`, `calculation::InputIssue`, `calculation::ChargeIssue`,
  `adapters::ImportIssue`, and `adapters::ResultDiagnostic` give their conformer, atom, and bond
  indices in one `core::IssueLocation location`, which also lists every atom and bond they concern
  in `atom_indices` and `bond_indices`.
- `features::PreparedMolecule` finds the geometry defects of its conformers once, in
  `geometry_defects()`, and copies share them. The defects include atoms that are too close, in
  `close_atoms`; `features::find_close_atoms()` finds those for one conformer without caching.
  `ConformerFeatures::coincident_atom_indices()` is removed.
- `core::Position` has a lexicographic three-way comparison.

### 3.0.0rc2

#### Python

- The execution mode and radius of a plan, an executed plan, and a rejection are in `execution`, an
  `Execution` value, matching the `execution` argument and the result JSON. Before, they were in
  `policy`, an `ExecutionPolicy` value.
- `Method.supports_fixed_point_sources` is renamed to `Method.supports_fixed_ions`.
- `Method.priority`, `ParameterSet.priority`, and `ParameterSetCatalog.for_method()` are removed. The
  catalogs list methods and parameter sets in the order automatic selection tries them; use
  `chargefw.methods["eem"].parameter_sets` for the sets of one method.
- The message of `NoExecutablePlanError` lists the first rejected methods and parameter sets with
  their reasons. Before, it was only `no_executable_plan`.
- Each fixed ion in `result.plan.fixed_ions.sources` has a `component_id`, the residue name it was
  selected by, such as `"MG"`, as in the result JSON.
- A calculation that runs out of memory raises `MemoryError`. Before, it raised
  `NumericalFailureError`.
- `Method.requires_bonds` is true for methods that use the bonds of the molecule.
- `Method.recommended` is true for methods whose charges are recommended for general use.
  `assess()` and `calculate()` take `method_scope="all"` to let automatic selection also consider the
  other methods; `result.requested.method_scope` records it.
- `chargefw.io.gemmi.from_structure()` accepts `bonds="hybrid"` and uses it by default: covalent links
  and disulfide bridges in `structure.connections` are added to the template bonds. Before, the only
  choices were `"templates"` and `"none"`, and the links were left out.

#### Command line

- `chargefw methods METHOD` reports `supports fixed ions` instead of `supports fixed point sources`.
  `chargefw methods METHOD` and `chargefw parameters ID` no longer print a priority; the listings are
  in the order automatic selection tries them.
- PEOE and MPEOE are named after partial equalization of *orbital* electronegativity, as in their
  publications.
- A calculation that runs out of memory reports `Not enough memory` and exits with status 1, without
  a result JSON. Before, it was reported as a numerical failure with status 4.
- `chargefw methods METHOD` reports whether the method requires bonds.
- `chargefw methods METHOD` reports whether the method is recommended. `calculate` and `applicability`
  take `--method-scope all` to let automatic selection also consider methods that are not
  recommended; the result JSON records it as `method_scope` in the requested provenance.

#### C++

- `calculate()` with an assessment no longer reports exceptions other than `std::logic_error` and
  `std::runtime_error` as `numerical_failure`; they reach the caller, for example `std::bad_alloc`.
- `MethodRequirements::bonds` marks methods that use the bonds of the molecule.
- `method_scope_from_string()` and `to_string()` convert `MethodScope` to and from its names.

#### Calculation

- Methods that use bonds (ABEEM, Charge2, DelRe, DENR, GDAC, KCM, MGC, MPEOE, PEOE, SQE, SQE+q0,
  SQE+qp, and TSEF) refuse a molecule of several atoms without bonds, as methods that use coordinates
  refuse a molecule without them. Before, some of them calculated meaningless charges for it.
- Automatic selection uses only methods whose charges are recommended for general use, so it no
  longer uses Charge2, DelRe, GDAC, SMP/QEq, VEEM, or the formal and dummy charges. Before, it fell
  back to them when no other method applied. They are still calculated when named.
- With permissive matching, an atom without bonds, such as the oxygen of a water molecule without its
  hydrogens, matches single-bonded parameters. Before, such an atom stopped the calculation even with
  permissive matching.

#### Input

- PDB and mmCIF records in which a residue of two or more atoms has an atom without bonds receive an
  `unbonded_residue_atoms` warning, for example a ligand in a PDB file without `CONECT` records.
- Residue templates bond an `H1` hydrogen to the nitrogen of an N-terminal proline, as Hydride names
  it. Before, reading such a structure with template bonds left the hydrogen unbonded.
- MOL and SDF V2000 records without `M  CHG` and `M  RAD` lines take formal charges from the charge field of the
  atom lines, as the format defines. Before, their atoms were read as neutral. A charge code outside
  0-7 is an error.

#### Documentation

- The new Python API reference lists every public class and attribute of the Python package.
- The method descriptions say which methods are not recommended for general use and why, and where
  ChargeFW differs from the publications, for example in the atom types of GDAC and KCM and the
  non-periodic EQeq.

### 3.0.0rc1

#### Command line

- A MOL2, mmCIF, PQR, or result JSON file that cannot be written no longer replaces an existing file of
  the same name. Before, a refused export, such as MOL2 for a molecule without coordinates, left the
  earlier file empty.
- A method option given as NaN or infinity, such as `--method-option denr.step=nan`, is refused as
  invalid input, exit status 2; it used to pass validation and fail later during the calculation.

#### Python

- `chargefw.io.rdkit.from_mol()` refuses an RDKit molecule with a 2D conformer. Embed 3D coordinates,
  remove the conformers, or mark a real flat geometry with `Set3D(True)`.
- `molecule.import_diagnostics` lists the warnings a reader attached to an imported molecule, such as
  `partial_charges_ignored`, as `chargefw.ImportDiagnostic` values with a code, a message, and a line.
  Before, they were only in the result JSON.
- `chargefw.calculate()` with an assessed plan refuses `parameter_matching`, `execution`,
  `cutoff_threshold`, and `cover_threshold` also when they are passed with their default value, such as
  `parameter_matching="strict"`. Before, such arguments were silently ignored and the plan's own
  settings were used.

#### C++

- Assessment refuses a custom method whose requirements list common parameters but no atom or bond
  parameters, with `std::invalid_argument`; declare such constants as method options. Common
  parameters belong in a parameter set only when they are fitted together with atom or bond
  parameters. Such a method used to pass assessment and then fail during execution.
  `MethodRequirements::requires_parameter_classification()` states which methods need atom or bond
  parameters.

#### Output

- Writing mmCIF for many molecules with the same name, such as an SDF library whose records share a
  title, no longer slows down with the cube of their number: 4,000 such molecules took about 100
  seconds and now take half a second. Result JSON and MOL2 output of large collections is also
  faster.
- ChargeFW writes PQR files for electrostatics with APBS: `--output-pqr` on the command line,
  `format="pqr"` in Python, and `PqrWriter` in C++. PQR is available for one structure read from PDB or
  mmCIF, with one conformer; each atom gets its element's van der Waals radius. ChargeFW refuses
  structures that APBS would misread, such as chain IDs other than a single letter.
  [File formats](https://chargefw.readthedocs.io/en/latest/reference/formats.html#pqr) lists the rules.
- mmCIF output contains the bonds used in the calculation, with their orders, so viewers that read
  mmCIF bond data can show double and triple bonds. Bonds within a residue are in `_chem_comp_bond`;
  bonds between residues, such as peptide links and disulfides, are in `_struct_conn`. ChargeFW refuses
  to write a structure whose bonds the residue definitions cannot express, for example two residues of
  the same name with different bonds.
  [File formats](https://chargefw.readthedocs.io/en/latest/reference/formats.html#mmcif) lists the rules.
- An mmCIF file can hold several charge sets for the same molecules: pass a list of results to
  `chargefw.io.dumps()`, `chargefw.io.write()`, `chargefw.io.gemmi.to_document()`, or
  `chargefw.io.gemmi.attach_charges()`. Two sets with the same method and parameter set are refused.
- A charge set calculated per conformer now has one ID in mmCIF output, covering all conformers. Before,
  each conformer got its own ID, so viewers listed the same method once per conformer.
- Attaching charges to an mmCIF document replaces its other charges dictionary declarations in
  `_audit_conform`, such as `mmcif_charges_v10.dic` from ACC or ChargeFW2 files, and repeated ones;
  declarations of other dictionaries stay.

#### Input

- An SDF or MOL record that declares a 2D drawing, such as a PubChem 2D download, is read without
  coordinates and gets a `coordinates_2d_ignored` warning. Methods that need 3D geometry no longer
  calculate charges from the flat drawing; the others still apply.
- MOL and SDF atoms no longer get their element symbol as a name, since molfiles have no atom names;
  `atom_names` in Python are empty for them. MOL2 and mmCIF output name these atoms from the element
  and atom position (`O1`, `O2`, `N4`), where they used to show bare symbols (`O`, `O`) or numbered
  repeats (`O`, `O2`).
- `chargefw.io.gemmi.from_structure()` uses template bonds by default and refuses `bonds="explicit"`
  and `"hybrid"`. Converting a `gemmi.Structure` loses the file's ligand bond orders, connection bond
  orders, and PDB `CONECT` records, so those bonds were lost or became single without a warning. With
  the new default, links that only the file records, such as disulfides and covalent ligands, are no
  longer bonded. For bonds from the file, use `from_document()` for mmCIF and `chargefw.io.read()` for
  PDB.
- PDB and mmCIF connections to an atom of a symmetry-related copy (`LINK`, `SSBOND`, or
  `_struct_conn` partners with different symmetry codes) are left out with a
  `symmetry_connections_ignored` warning. Before, they became bonds to the atom of the same name in the
  file, which is a different atom.
- An mmCIF block that contains partial charges, for example from an earlier calculation, gives its
  record a `partial_charges_ignored` warning, as MOL2 input does; the charges are not used.
- mmCIF bonds between residues keep their order from `_struct_conn.pdbx_value_order`; they were all
  read as single bonds before. A connection without an order is still a single bond.

### 3.0.0a4

#### Installation

- ChargeFW requires Python 3.12 or newer. One wheel serves all of them, including Python versions
  released later. On Python 3.11, pip installs an older ChargeFW release instead;
  [Installation](https://chargefw.readthedocs.io/en/latest/getting-started/installation.html) shows how
  to create an environment with a newer Python.
- `pip install chargefw` also installs the `chargefw` command line.
- Wheels are also published for Linux aarch64, macOS 14 or newer on Apple silicon, and Windows
  x86-64. On Windows, the command line accepts file names that the system code page cannot represent.

#### Distribution

- The wheel contains the ChargeFW library and Gemmi in its Python extension and command line, and no
  longer ships separate `libchargefw_core.so` and `libgemmi_cpp.so` libraries. It is about a quarter
  smaller.

#### Building from source

- The new CMake option `CHARGEFW_BUILD_LIBRARY`, on by default, builds and installs the shared C++
  library. The library and its installation are supported on Linux; on Windows, CMake stops with an
  explanation and points to the Python package instead.
- ChargeFW compiles with Clang and libc++ 19 or newer, including Apple Clang from Xcode 16.4 and 26.1.

#### Python

- The type aliases in `chargefw.core` and `chargefw.io`, such as `chargefw.io.InputFormat`, are defined
  with the `type` statement. Type checkers treat them as before; at run time, `typing.get_args()` needs
  the alias's `__value__`, for example `typing.get_args(chargefw.io.InputFormat.__value__)`.

### 3.0.0a3

#### Documentation

- The [documentation](https://chargefw.readthedocs.io/en/latest/) is reorganized and rewritten for
  users of the Python package, the command line, and the C++ library. Pages have new addresses, so
  links to the old ones, such as `PYTHON.html`, `FORMATS.html`, or `recipes/`, no longer work.
  [Differences from ChargeFW2](https://chargefw.readthedocs.io/en/latest/about/chargefw2.html)
  describes what changed from ChargeFW2 and how to move over.

#### Output

- MOL2 and mmCIF files write charges rounded to four decimal places, and coordinates with the digits
  they were read with instead of 17 significant digits. Result JSON keeps full precision.

#### Listings

- `chargefw methods`, `chargefw.methods`, and the C++ method registry list methods by descending
  priority, then ID. `chargefw parameters`, `chargefw.parameter_sets`, and the C++ parameter loaders list
  parameter sets grouped by method, each method's sets by descending priority, then ID.

### 3.0.0a2

#### Building from source

- Native builds use installed CLI11, nlohmann/json, Eigen, nanoflann, oneTBB, and Gemmi packages by
  default and report any that are missing at configure time. The supported versions are listed in
  the [native library documentation](https://chargefw.readthedocs.io/en/latest/interfaces/cpp.html#build-and-link);
  the minimums match Debian 13. Set
  `CHARGEFW_FETCH_DEPENDENCIES=ON` to download pinned sources instead, as the developer presets and
  the Python package build do. This option replaces `CHARGEFW_USE_SYSTEM_DEPENDENCIES`; builds that
  set it to `ON` should use `CHARGEFW_FETCH_DEPENDENCIES=OFF`.
- Downloaded dependency sources are verified against SHA-256 checksums.
- The installed CMake package requires the exact nlohmann/json and Gemmi versions that ChargeFW was
  built against.
- Installed binaries find dependency libraries installed outside the system library paths without
  `LD_LIBRARY_PATH`.

#### Distribution

- Wheels, native installations, and the CLI image include the license notices of the third-party
  code in their binaries.

#### Fixes

- mmCIF input accepts tokens longer than 4096 bytes; entries with long values, such as 3J3Q, previously
  failed to parse.
- Cover execution of SQE+qp conserves the total charge on very large structures; on 3J3Q (2.44 million
  atoms), rounding had left the total at 1.1e-4 instead of 0.

### 3.0.0a1

First pre-release of the ChargeFW 3 series, the successor of
[ChargeFW2](https://github.com/sb-ncbr/ChargeFW2). The version number continues ChargeFW2's lineage;
the `chargefw` Python package moves from 0.1.7 directly to 3.0.0a1.

#### Changes since 0.1.7

- Parameter-set keys are exposed in Python (`atom_keys` and `bond_keys`) and in the CLI
  (`chargefw parameters ID --keys`).
- The bundled periodic table is exposed in Python as `chargefw.periodic_table`.
- New recipes calculate dipole moments and attach charges to RDKit SDF records; the file calculation
  recipe can export CSV, and the fixed-ion and molecule-inspection recipes report more detail.
- The repository moved to [sb-ncbr/chargefw](https://github.com/sb-ncbr/chargefw).
