# Changelog

## Unreleased

### Installation

- On x86-64, the wheels, the Docker image, and source builds target x86-64-v3 (AVX2 and FMA: Intel
  Haswell, AMD Excavator, or newer). Full calculations on large structures run about twice as fast
  with the wheels. On an older CPU, `import chargefw` raises an `ImportError`; build from source with
  `CHARGEFW_ENABLE_NATIVE_OPTIMIZATIONS=ON`, which now defaults to `OFF`.

### Python

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

### C++

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
- `adapters::to_string()` names a `MolecularSourceFormat` and a `SourceConnectivity`, and
  `molecular_source_format_from_string()` reads a format name.
- `calculation::calculate(assessment, plan, candidate)` runs a plan with other values for its parameter
  set, without checking the molecules again.
- `parameters::save_parameter_set_json()` writes a parameter set in the JSON format that the loading
  functions read.
- `MethodRegistry` cannot be constructed; `method_registry()` is the catalog of built-in methods.
  `find_applicable_methods()`, `ApplicabilityRequest`, `CalculationRequest`, `CalculationResult`, and
  the `calculate()` overload that takes a `CalculationRequest` are removed. Calculate charges with
  `assess()` and `calculate()`.

### Calculation

- SQE, SQE+q0, and SQE+qp calculations fail with `numerical_failure` when the split-charge system
  cannot be factorized as positive definite, as when the energy has no unique minimum at the
  molecule's geometry. Earlier versions returned the solution without warning, with charges that could
  be far too large. In cutoff and cover execution, each fragment is checked. Solving the system takes
  about half the time.

## 3.0.0rc3

### Input

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

### Calculation

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

### Output

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

### Python

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

### Command line

- An unexpected error during assessment is reported as an internal error, exit status 1. Before, it
  was reported as a numerical failure, exit status 4, with result JSON.
- Each warning of the executed plan, such as a full calculation above the size limit, is printed once,
  before the calculation. Before, it was printed once per molecule.
- A warning about a molecule names it, such as `Warning: molecule 1 ('1CRN'): …`. Several molecules
  with the same warning give one line with their number and the first of them. Before, the line did
  not say which molecule it concerned.

### C++

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

## 3.0.0rc2

### Python

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

### Command line

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

### C++

- `calculate()` with an assessment no longer reports exceptions other than `std::logic_error` and
  `std::runtime_error` as `numerical_failure`; they reach the caller, for example `std::bad_alloc`.
- `MethodRequirements::bonds` marks methods that use the bonds of the molecule.
- `method_scope_from_string()` and `to_string()` convert `MethodScope` to and from its names.

### Calculation

- Methods that use bonds (ABEEM, Charge2, DelRe, DENR, GDAC, KCM, MGC, MPEOE, PEOE, SQE, SQE+q0,
  SQE+qp, and TSEF) refuse a molecule of several atoms without bonds, as methods that use coordinates
  refuse a molecule without them. Before, some of them calculated meaningless charges for it.
- Automatic selection uses only methods whose charges are recommended for general use, so it no
  longer uses Charge2, DelRe, GDAC, SMP/QEq, VEEM, or the formal and dummy charges. Before, it fell
  back to them when no other method applied. They are still calculated when named.
- With permissive matching, an atom without bonds, such as the oxygen of a water molecule without its
  hydrogens, matches single-bonded parameters. Before, such an atom stopped the calculation even with
  permissive matching.

### Input

- PDB and mmCIF records in which a residue of two or more atoms has an atom without bonds receive an
  `unbonded_residue_atoms` warning, for example a ligand in a PDB file without `CONECT` records.
- Residue templates bond an `H1` hydrogen to the nitrogen of an N-terminal proline, as Hydride names
  it. Before, reading such a structure with template bonds left the hydrogen unbonded.
- MOL and SDF V2000 records without `M  CHG` and `M  RAD` lines take formal charges from the charge field of the
  atom lines, as the format defines. Before, their atoms were read as neutral. A charge code outside
  0-7 is an error.

### Documentation

- The new Python API reference lists every public class and attribute of the Python package.
- The method descriptions say which methods are not recommended for general use and why, and where
  ChargeFW differs from the publications, for example in the atom types of GDAC and KCM and the
  non-periodic EQeq.

## 3.0.0rc1

### Command line

- A MOL2, mmCIF, PQR, or result JSON file that cannot be written no longer replaces an existing file of
  the same name. Before, a refused export, such as MOL2 for a molecule without coordinates, left the
  earlier file empty.
- A method option given as NaN or infinity, such as `--method-option denr.step=nan`, is refused as
  invalid input, exit status 2; it used to pass validation and fail later during the calculation.

### Python

- `chargefw.io.rdkit.from_mol()` refuses an RDKit molecule with a 2D conformer. Embed 3D coordinates,
  remove the conformers, or mark a real flat geometry with `Set3D(True)`.
- `molecule.import_diagnostics` lists the warnings a reader attached to an imported molecule, such as
  `partial_charges_ignored`, as `chargefw.ImportDiagnostic` values with a code, a message, and a line.
  Before, they were only in the result JSON.
- `chargefw.calculate()` with an assessed plan refuses `parameter_matching`, `execution`,
  `cutoff_threshold`, and `cover_threshold` also when they are passed with their default value, such as
  `parameter_matching="strict"`. Before, such arguments were silently ignored and the plan's own
  settings were used.

### C++

- Assessment refuses a custom method whose requirements list common parameters but no atom or bond
  parameters, with `std::invalid_argument`; declare such constants as method options. Common
  parameters belong in a parameter set only when they are fitted together with atom or bond
  parameters. Such a method used to pass assessment and then fail during execution.
  `MethodRequirements::requires_parameter_classification()` states which methods need atom or bond
  parameters.

### Output

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

### Input

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

## 3.0.0a4

### Installation

- ChargeFW requires Python 3.12 or newer. One wheel serves all of them, including Python versions
  released later. On Python 3.11, pip installs an older ChargeFW release instead;
  [Installation](https://chargefw.readthedocs.io/en/latest/getting-started/installation.html) shows how
  to create an environment with a newer Python.
- `pip install chargefw` also installs the `chargefw` command line.
- Wheels are also published for Linux aarch64, macOS 14 or newer on Apple silicon, and Windows
  x86-64. On Windows, the command line accepts file names that the system code page cannot represent.

### Distribution

- The wheel contains the ChargeFW library and Gemmi in its Python extension and command line, and no
  longer ships separate `libchargefw_core.so` and `libgemmi_cpp.so` libraries. It is about a quarter
  smaller.

### Building from source

- The new CMake option `CHARGEFW_BUILD_LIBRARY`, on by default, builds and installs the shared C++
  library. The library and its installation are supported on Linux; on Windows, CMake stops with an
  explanation and points to the Python package instead.
- ChargeFW compiles with Clang and libc++ 19 or newer, including Apple Clang from Xcode 16.4 and 26.1.

### Python

- The type aliases in `chargefw.core` and `chargefw.io`, such as `chargefw.io.InputFormat`, are defined
  with the `type` statement. Type checkers treat them as before; at run time, `typing.get_args()` needs
  the alias's `__value__`, for example `typing.get_args(chargefw.io.InputFormat.__value__)`.

## 3.0.0a3

### Documentation

- The [documentation](https://chargefw.readthedocs.io/en/latest/) is reorganized and rewritten for
  users of the Python package, the command line, and the C++ library. Pages have new addresses, so
  links to the old ones, such as `PYTHON.html`, `FORMATS.html`, or `recipes/`, no longer work.
  [Differences from ChargeFW2](https://chargefw.readthedocs.io/en/latest/about/chargefw2.html)
  describes what changed from ChargeFW2 and how to move over.

### Output

- MOL2 and mmCIF files write charges rounded to four decimal places, and coordinates with the digits
  they were read with instead of 17 significant digits. Result JSON keeps full precision.

### Listings

- `chargefw methods`, `chargefw.methods`, and the C++ method registry list methods by descending
  priority, then ID. `chargefw parameters`, `chargefw.parameter_sets`, and the C++ parameter loaders list
  parameter sets grouped by method, each method's sets by descending priority, then ID.

## 3.0.0a2

### Building from source

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

### Distribution

- Wheels, native installations, and the CLI image include the license notices of the third-party
  code in their binaries.

### Fixes

- mmCIF input accepts tokens longer than 4096 bytes; entries with long values, such as 3J3Q, previously
  failed to parse.
- Cover execution of SQE+qp conserves the total charge on very large structures; on 3J3Q (2.44 million
  atoms), rounding had left the total at 1.1e-4 instead of 0.

## 3.0.0a1

First pre-release of the ChargeFW 3 series, the successor of
[ChargeFW2](https://github.com/sb-ncbr/ChargeFW2). The version number continues ChargeFW2's lineage;
the `chargefw` Python package moves from 0.1.7 directly to 3.0.0a1. Install it with
`pip install --pre chargefw`.

ChargeFW 3 is a new implementation. Its highlights are:

- one native engine behind a C++23 library, a command-line application, and a Python package for
  CPython 3.11 to 3.14 on Linux x86-64;
- 22 methods and the bundled parameter sets, with deterministic method and parameter-set selection and
  structured reasons for rejected candidates;
- explicit molecular input: topology, protonation, formal charges, and coordinates are never changed
  silently, and missing prerequisites are reported as failures;
- full calculations and explicit cutoff and cover approximations, recorded in result provenance;
- SDF, MOL2, PDB, mmCIF, and molecule JSON input; result JSON, MOL2, and mmCIF output; source atom,
  molecule, and conformer order preserved in every result;
- optional Gemmi and RDKit integration in Python.

Documentation is published at [chargefw.readthedocs.io](https://chargefw.readthedocs.io/).

### Changes since 0.1.7

- Parameter-set keys are exposed in Python (`atom_keys` and `bond_keys`) and in the CLI
  (`chargefw parameters ID --keys`).
- The bundled periodic table is exposed in Python as `chargefw.periodic_table`.
- New recipes calculate dipole moments and attach charges to RDKit SDF records; the file calculation
  recipe can export CSV, and the fixed-ion and molecule-inspection recipes report more detail.
- The repository moved to [sb-ncbr/chargefw](https://github.com/sb-ncbr/chargefw).
