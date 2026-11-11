# ChargeFW code review

A read-only review of the native library, the command-line application, and the Python package, made
on 2026-11-10 at commit `cff64e0` (version 3.0.0). It looked for bugs, lifetime problems, public API
inconsistencies, dead code, and duplicated code that justify small, targeted changes. Style, naming
taste, micro-optimizations, and large redesigns were out of scope.

## How to read this file

Each finding has an ID, so that it can be picked up on its own. Before fixing one, reassess it: the
code may have changed, and most findings were established by reading only.

- **Nothing was built or run for this review.** No finding was reproduced at runtime, except where a
  finding says so. Reproduction snippets are suggestions that have not been executed.
- **Line numbers refer to commit `cff64e0`.** If they no longer match, search for the quoted code.
- **Status** says how far a finding was checked:
  - *re-checked*: found by one reader and then confirmed against the source by a second read;
  - *reported*: found by one reader, not independently confirmed;
  - *plausible*: depends on something that was not checked, which the finding names.
- **Check** gives a command or a short procedure that confirms or refutes the finding.
- **Fix** is the smallest change that was considered sufficient. It is a proposal, not a decision.
- **Resolution** appears where a code change is not the obvious answer. It says whether to fix the
  code, document the behaviour, or leave it.

### Fix in code or document

Not every finding should be fixed in code. The rule used here:

- **Fix in code** when a common workflow reaches the problem: the command line, the Python package,
  or the C++ path `assess()` → `calculate()` → writers, with files that real programs produce.
- **Document** when only an unusual use reaches it, such as calling a lower-level C++ function
  directly or using an object in a way no workflow does, and a sentence tells the reader what to do
  instead. A notice costs less than a new check, a new failure path, and its test.
- **Leave** when it needs deliberate misuse, or when nothing a reader would do leads to it.

A fix that is as small as the notice and that the compiler enforces still wins over the notice.

Where a notice goes follows the document ownership in AGENTS.md: a precondition of a C++ function
goes into the comment at its declaration; behaviour a C++ user must know goes into
`docs/interfaces/cpp.md`; file-format behaviour goes into `docs/reference/formats.md`. User
documentation mentions a limitation only where the reader needs it to avoid a wrong result or a
failed run, together with what to do instead.

Findings where this rule changes the proposal:

| Finding | Resolution | Where |
| --- | --- | --- |
| L1 | Already documented; code change optional | `docs/interfaces/cpp.md:5` |
| L2 | Document | comment in `prepared_molecule.h` |
| B6 | Fix or document, after checking real files | `docs/reference/formats.md` |
| B8 | Document | `molecule_record.h`, result schema, `docs/reference/formats.md` |
| B10 | Leave, or delete the `reserve` calls | — |
| A6 | Delete the unused constructor only | — |
| A7 | Document the stream check | comment in `json_output.h` |
| A10 | Leave | — |
| A12 (enum members) | Leave | — |
| X9 (PDB serials, trimmed V2000 lines) | Document if confirmed | `docs/reference/formats.md` |

Every other finding in sections B and A is reached by a common workflow, or is a consistency or
cleanup change with no behaviour to document, and stays a code change.

Follow [AGENTS.md](AGENTS.md) when acting on a finding. In particular, the findings marked
*design-sensitive* change public API or layering and need a sketch and the maintainer's review before
implementation. Remove a finding from this file when it is fixed or rejected; do not keep history
here.

The greps below are meant to be run from the repository root. "Repo-wide" means
`include src apps python tests`.

## Summary

No memory-safety bug and no data race was found in the library. The most valuable fixes are the
behaviour bugs in section B, all of which give a wrong result, a failed run, or a misleading message
for realistic input.

| Section | Content | Findings |
| --- | --- | --- |
| B | Behaviour bugs | B2–B11 |
| L | Lifetime and misuse traps | L1–L2 |
| A | Public API design | A3–A12 |
| D | Dead code | D2–D13 |
| C | Duplicated code | C1–C15 |
| X | Checked and found clean, or rejected | X1–X9 |

---

## B. Behaviour bugs

### B2. MOL2 silently drops atoms or bonds when the declared count is too small

- **Status:** re-checked.
- **Where:** `src/adapters/native/mol2_input.cpp:85` (`read_atoms`), `:144` (`read_bonds`), `:59-67`
  (`next_section`), `:234-246` (section loop).
- **Problem:** `read_atoms` and `read_bonds` read exactly the declared number of lines. `next_section`
  then discards every following line that does not start with `@`, so surplus atom or bond lines are
  skipped as if they were padding. A count that is too large is an error; a count that is too small
  is not.
- **Failure:** a counts line `8 5` followed by seven bond lines imports with two bonds missing and no
  warning. `docs/reference/formats.md` states that wrong atom or bond counts are errors.
- **Check:** copy a small MOL2 fixture, lower the bond count on the counts line by one, and read it;
  the record is returned with one bond fewer instead of an error.
- **Fix:** after `read_atoms` and `read_bonds`, require the next non-blank, non-comment line to start
  with `@` or be the end of input; otherwise throw `std::runtime_error`, as the reader does for other
  malformed input.
- **Tests:** one case in `tests/cpp/adapters/` next to the existing MOL2 count checks.

### B3. Writing result JSON throws on non-UTF-8 text taken from input files

- **Status:** re-checked for the code path; the CLI message follows from the catch order.
- **Where:** `src/adapters/native/json_output.cpp:321`; catch order in `apps/chargefw/main.cpp:270-285`.
- **Problem:** `document.dump(2)` uses nlohmann's strict error handler and throws
  `nlohmann::json::type_error` (id 316) for invalid UTF-8. Record IDs, molecule names, the source
  path, and diagnostic messages are copied from input files without re-encoding.
- **Failure:** an SDF whose title line is Latin-1 (for example `Benz\xE8ne`) calculates successfully.
  Writing the result then throws an exception that is neither `std::runtime_error` nor
  `std::invalid_argument`, so the CLI prints "Internal error" and exits with status 1.
- **Check:** write an SDF fixture with byte `0xE8` in the title and run the installed CLI through the
  registered CTest CLI machinery (not from the build tree; see AGENTS.md).
- **Fix:** `document.dump(2, ' ', false, nlohmann::json::error_handler_t::replace)`.
- **Tests:** one case in `tests/cpp/adapters/test_json_output.cpp` with an invalid byte in a record ID.

### B4. Ctrl+C is ignored outside the calculation loop

- **Status:** re-checked.
- **Where:** `apps/chargefw/main.cpp:27-31` (flag and handler), `:54-77` (`ScopedInterruptHandler`),
  `:218` (its scope), `:261` (`write_calculation_outputs`).
- **Problem:** the handler is installed at `:218` and lives until `run()` returns, so it is still
  installed while outputs are written, where nothing reads the flag. The handler only sets the flag;
  it never restores the default action, so a second Ctrl+C does nothing either.
- **Failure:** on POSIX the process cannot be stopped with Ctrl+C during `assess()`, during one long
  full solve (a single target is not interruptible), or while large outputs are written.
  `docs/interfaces/cli.md:72-76` says writing does not check for cancellation, but not that the
  signal is swallowed.
- **Check:** read `main.cpp:210-265`; `interrupt_handler` is declared before the `execution` lambda
  and destroyed at the end of the function.
- **Fix:** move `interrupt_handler` and `interruptible_observer` inside the `execution` lambda. Make a
  second signal fatal: in `request_interruption`, if `test_and_set` returns true, call
  `std::signal(SIGINT, SIG_DFL)` and `std::raise(SIGINT)`. Update `docs/interfaces/cli.md` if the
  described behaviour changes.

### B5. Fixed-ion errors print zero-based indices

- **Status:** re-checked (lines 24-41 and 100-129 were read).
- **Where:** `src/calculation/fixed_charge_partition.cpp:27-28, 33-35, 39-40, 50-53, 76-78, 82-84,
  89-90, 94-95, 105-107, 126-129`.
- **Problem:** these messages print raw indices. Every other diagnostic is one-based through
  `core::detail::molecule_description` and `atom_description`
  (`src/core/diagnostic_description.h:31, 52`). The same feature mixes both conventions:
  `remap_fixed_ion_issue` in `src/calculation/assessment.cpp:214, 226` is one-based.

  ```cpp
  "fixed charge source at molecule " + std::to_string(source.molecule_index) +
  ", atom " + std::to_string(source.atom_index) + " has non-finite charge"
  ```
- **Failure:** `chargefw calculate --fixed-ions CA` on a structure with a bonded ion names
  "molecule 0, atom 1234"; a reader who treats it like any other message looks at the wrong atom.
- **Check:** `grep -n "std::to_string(.*_index)" src/calculation/fixed_charge_partition.cpp`.
- **Open question:** the zero-based text may have been meant to mirror `FixedAtomCharge` indices for
  native callers. Decide before changing.
- **Fix:** build the messages with `molecule_description` / `atom_description`.
- **Tests:** `tests/cpp/calculation/test_fixed_ion_planning.cpp:692, 712` pin the current wording and
  need updating; per AGENTS.md they should not assert complete wording at all.

### B6. V3000 reader rejects valid records

- **Status:** reported. That RDKit omits an empty bond block is from memory, not checked.
- **Where:** `src/adapters/native/mol_input.cpp:321-324` and `:368-372`.
- **Problem:** after the atom block the reader requires exactly `END ATOM`, `BEGIN BOND`, and after
  the bond block exactly `END BOND`, `END CTAB`, `M  END`. In the CTfile format the bond block is
  optional when there are no bonds, and SGROUP, COLLECTION, and 3D blocks may follow it.
- **Failure:** a single-atom V3000 record without a bond block, or a molecule with
  `M  V30 BEGIN COLLECTION` (stereo groups), is rejected with "expected V3000 bond block".
- **Check:** write a V3000 record with one atom and no `BEGIN BOND`/`END BOND` pair and read it. The
  tests cover only an empty bond block.
- **Fix:** after `END ATOM`, loop over V30 lines: parse `BEGIN BOND`, skip any other
  `BEGIN x` … `END x` block, stop at `END CTAB`.
- **Resolution:** first check what real writers produce (RDKit and Open Babel V3000 output for a
  single atom and for a molecule with stereo groups). If they produce such records, fix the code:
  SDF from those programs is a common workflow. If not, state in `docs/reference/formats.md`, next
  to the V3000 description, that a record needs a bond block and that blocks after it are rejected.

### B8. `ImportIssue::line` restarts at each SDF record

- **Status:** reported.
- **Where:** `src/adapters/native/mol_input.cpp:423` (`std::size_t line = 0;` per `parse_mol` call),
  used at `:179`, `:356`, `:398`.
- **Problem:** `include/chargefw/adapters/molecule_record.h:203` documents "One-based line in the
  source file" and `schemas/result-1.0.schema.json:227` says "One-based source line number". For SDF
  records after the first, the number counts from the start of the record.
- **Check:** read an SDF with an ignored V2000 property in its second record and inspect the issue's
  line. Tests assert lines only for first records.
- **Fix:** let `SdfReader` keep a running line count (an internal `parse_mol` overload taking
  `std::size_t& line`, and count lines in `consume_to_sdf_delimiter`). The alternative is to change
  both descriptions to "line within the record".
- **Resolution:** document. The number is useful either way, and the record index is reported with
  it. Change the comment in `molecule_record.h`, the description in the result schema, and the
  sentence about line numbers in `docs/reference/formats.md` to say that the line is counted within
  the record.

### B9. MOL2 section markers are matched exactly, and repeated sections overwrite

- **Status:** reported. How common trailing whitespace is in real files is unknown.
- **Where:** `src/adapters/native/mol2_input.cpp:234-246`, `:333`.
- **Problem:** only `\r` is stripped before a marker is compared, and a second ATOM or BOND section
  replaces the first.
- **Failure:** a second record whose marker is `@<TRIPOS>MOLECULE ` with a trailing space is treated
  as an unknown section. Its ATOM and BOND sections then overwrite those of the first record when
  the counts match, as in multi-conformer files, and the second record is lost. A trailing space on
  the first marker yields no records.
- **Fix:** compare `common::trim(*section)`, and throw on a repeated ATOM or BOND section within one
  record.

### B10. Declared counts are reserved before any line is read

- **Status:** plausible; the outcome depends on the allocator and overcommit settings.
- **Where:** `src/adapters/native/mol_input.cpp:275-278` (V3000),
  `src/adapters/native/mol2_input.cpp:81-83`, `:142`.
- **Problem:** the declared count is checked only for being positive before vectors are reserved.
- **Failure:** a 200-byte file with `M  V30 COUNTS 2000000000 0`, or a MOL2 counts line
  `2000000000 0`, requests tens of gigabytes and ends in `std::bad_alloc`, which the CLI reports as
  "Not enough memory" with status 1 instead of invalid input. V2000 is bounded at 999 by its field
  width.
- **Fix:** drop these `reserve` calls or cap them, for example `std::min(count, 1 << 16)`.
- **Resolution:** leave, or delete the `reserve` calls if the file is being edited anyway. Only a
  deliberately malformed file reaches this, and the run still fails with a message. Nothing to
  document.

### B11. The JSON writer's component ID differs from the one fixed ions match

- **Status:** reported.
- **Where:** `src/adapters/native/json_output.cpp:104-117` (`component_id`) versus
  `detail::canonical_component_name` in `src/adapters/component.h:49-55`.
- **Problem:** the JSON copy tests only for a non-empty label. The shared helper also treats `.` and
  `?` as unknown, and the mmCIF reader stores those tokens verbatim
  (`src/adapters/gemmi/mmcif_labels.h:28-30`).
- **Failure:** for a file with `label_comp_id` `.` and `auth_comp_id` `NA`, fixed-ion resolution
  matches `NA` but the result JSON reports `"fixed_ions[].component_id": "."`.
- **Fix:** call `detail::canonical_component_name` from `json_output.cpp` and delete the local copy.

---

## L. Lifetime and misuse traps

### L1. `ExecutionPlan` is copyable but points into its assessment

- **Status:** re-checked.
- **Where:** `include/chargefw/calculation/assessment.h:43-77`.
- **Problem:** the class holds `std::shared_ptr<const PlanIdentity> identity_` and a raw
  `const methods::ApplicableMethod* candidate_` into `AssessmentResult::applicability_`. Copying is
  not disabled, and the shared pointer suggests that a copy stands on its own. Only the class comment
  says otherwise.
- **Failure:** `auto plan = *calculation::assess(request).default_plan(); plan.method().id();`
  compiles and reads freed memory.
- **Check:** `grep -rn "ExecutionPlan" src apps python/src tests/cpp | grep -v "const.*ExecutionPlan&\|ExecutionPlan\*"`
  to confirm that nothing copies a plan by value. The first reader found no copy; the Python binding
  returns `const ExecutionPlan&` (`python/src/calculation.h:50`), and the only construction is
  `src/calculation/assessment.cpp:394`.
- **Fix:** delete the copy constructor and copy assignment; default the move operations (`plans_`
  needs only move).
- **Resolution:** already documented. `docs/interfaces/cpp.md:5` names the assessment as the one
  lifetime to manage and says that its plans refer to it, and the class comment says the same. The
  documented workflow keeps the assessment alive. Deleting the copy operations is optional
  hardening: it is four lines, no caller is affected, and the compiler then enforces what the
  documentation asks for.
- **Validation:** `gcc-debug` and `clang-debug` builds, since this is overload- and
  conversion-sensitive; the Python extension must still compile.

### L2. A moved-from `PreparedMolecule` indexes out of bounds

- **Status:** reported. No code was found that uses a moved-from object.
- **Where:** `include/chargefw/features/prepared_molecule.h:58-71`,
  `src/features/topology_features.cpp:39-60`.
- **Problem:** the hand-written move operations share the geometry cache "so that geometry queries
  remain valid on the moved-from object", but move `topology_`. The moved-from `TopologyFeatures`
  keeps its molecule pointer while its vectors are emptied, so `degree`, `neighbor_indices`,
  `incident_bond_indices`, and `bond_index_between` pass `validate_atom_index` (which checks the
  molecule's atom count) and then index an empty vector.
- **Fix:** validate against `neighbor_indices_.size()`. The alternative is to copy `topology_` in the
  move operations, as is done for the cache.
- **Resolution:** document. No workflow uses a moved-from prepared molecule. Correct the comment
  above the move constructor in `prepared_molecule.h`: only `geometry_defects()` stays valid on the
  moved-from object; its topology must not be queried. This is not user documentation.

---

## A. Public API design

### A3. `RequestedCalculationProvenance` uses strings where enums exist (design-sensitive)

- **Status:** re-checked.
- **Where:** `include/chargefw/adapters/charge_result.h` (`method_scope = "recommended"`,
  `execution_kind = "auto"`); writers of the fields in `apps/chargefw/output.cpp:171, 176` and
  `python/src/calculation.cpp:408, 413`; reader in `src/adapters/native/json_output.cpp:253, 260`.
- **Problem:** `calculation::MethodScope` and `calculation::ExecutionSelectionKind` exist. The CLI
  converts enum to string; the Python binding passes the caller's string into provenance and parses
  the same string into the enum separately (`python/src/calculation.cpp:441`), so nothing ties the
  recorded value to the one that was used.
- **Fix:** store the enums and convert once in the JSON writer.

### A4. `methods` depends on a `calculation` header (design-sensitive)

- **Status:** re-checked.
- **Where:** `include/chargefw/methods/method_applicability.h:3` includes
  `chargefw/calculation/execution_policy.h` for `calculation::ExecutionMode`;
  `src/methods/find_applicable_methods.h:25` and `src/methods/method_applicability.cpp:39-40` use
  `calculation::ResourcePolicy`.
- **Problem:** the lower layer includes a header of the layer built on top of it.
- **Check:** `grep -rn "calculation::" src/methods include/chargefw/methods`.
- **Fix:** not proposed here. The options are to move `ExecutionMode` and `ResourcePolicy` down, or to
  move the execution-assessment part of applicability up into `calculation`. Sketch first.

### A5. `PrerequisiteResult` differs from its sibling result types

- **Status:** reported.
- **Where:** `include/chargefw/methods/method_prerequisites.h:39-49` (class with private `issues_`,
  `issues()` returning a span, `add()`), `include/chargefw/methods/parameter_prerequisites.h:34-41`
  (struct with a public `issues` vector). Both have `explicit operator bool`.
- **Cost:** `copy_issues` in `src/methods/method_applicability.cpp:55-58` exists only to bridge the
  two (`:218` versus `:260`).
- **Fix:** make `PrerequisiteResult` a struct with a public `std::vector<PrerequisiteIssue> issues`,
  or add `take_issues() &&`. Either removes `copy_issues`.

### A6. `ClassificationResult` can be built in a meaningless success state

- **Status:** reported.
- **Where:** `include/chargefw/parameters/classification/classification_result.h:29-32`,
  `src/parameters/classification/classification_result.cpp:11-16`.
- **Problem:** the default constructor, and the issues constructor given an empty vector, both yield
  a truthy result with an empty classification. `classify_parameters` relies on "falsy implies
  non-empty issues" when it calls `result.issues().front()`
  (`src/parameters/classification/parameter_classifier.cpp:305`).
- **Check:** `grep -rn "ClassificationResult{}\|ClassificationResult result;" include src apps python tests`;
  the first reader found no user of the default constructor.
- **Fix:** delete the default constructor; throw `std::invalid_argument` from the issues constructor
  when the vector is empty.
- **Resolution:** delete the unused default constructor, which is dead code. Do not add the throw:
  the only code that constructs the type is the classifier itself, so the new failure path would
  guard nothing.

### A7. `JsonWriter::write` does not check stream failure

- **Status:** reported.
- **Where:** `include/chargefw/adapters/native/json_output.h:16`.
- **Problem:** It is the only writer that does not throw when the stream fails; the others end with
  `if (!*output_) throw`. The CLI checks the stream itself, so only direct C++ users are affected.
- **Fix:** add the same stream check.
- **Resolution:** either add the check or say in the comment at the declaration that the caller checks
  the stream after `write()`; the example in `docs/interfaces/cpp.md:162` should then show that check.

### A8. Python hard-codes the resource thresholds

- **Status:** re-checked.
- **Where:** `python/chargefw/calculation.py:298-299, 393-394`,
  `python/chargefw/_calculation_options.py:152-153`.
- **Problem:** `20_000` and `80_000` are literals. The native library owns them as
  `calculation::default_cutoff_atom_threshold` and `default_cover_atom_threshold`
  (`include/chargefw/calculation/execution_policy.h:10-11`), and the CLI uses the constants
  (`apps/chargefw/input.cpp:257, 262`). A native change would leave Python on a different execution
  policy.
- **Check:** `grep -rn "20_000\|80_000" python/chargefw`.
- **Fix:** export the two constants from `bind_calculation` (as `_fixed_ion_names` is exported) and
  use them at the three default sites; update the `.pyi` stub.

### A10. "Read-only" fit arrays can be made writable again

- **Status:** reported; the NumPy behaviour was reproduced with plain NumPy.
- **Where:** `python/chargefw/parameterization/_evaluator.py:76-77` (`Evaluator.values`),
  `_fit.py:332, 447-448, 625-628` (`best_values`, `values`, `start`, `bounds`),
  `_diagnostics.py:162` (`singular_values`).
- **Problem:** these call `setflags(write=False)` on arrays that own their data, which a caller can
  reverse. `evaluator.values.setflags(write=True); evaluator.values[:] = 0` changes
  `Evaluator._values`, which later fits read (`_fit.py:567`).
- **Fix:** use `immutable_array()` from `python/chargefw/_arrays.py:136`, as the rest of the package
  does.
- **Resolution:** leave. Reaching this takes a deliberate `setflags(write=True)` on an array the
  package marked read-only. Switch to `immutable_array()` only if these lines are edited for another
  reason. Nothing to document.

### A11. Smaller Python API items

- **Status:** reported.
- `python/chargefw/calculation.py:132-148`: the `__module__` rewrite loop omits `InputIssue` and
  `ChargeIssue`, so they report `chargefw._calculation_values`. Add the two names.
- `python/chargefw/calculation.py:198`: `write_parameter_set(parameter_set, path)` has the reverse
  argument order of `io.write(path, result)`. Changing it would break the public API; this is only
  noted, pending a separate decision.
- `python/chargefw/io/gemmi.py:46, 145`: `StructureBondStrategy` is used in the `from_structure`
  signature but missing from `__all__`.
- `python/chargefw/_catalog.py:42-44`: `5 in chargefw.methods` raises `TypeError`, because
  `Mapping.__contains__` catches only `KeyError`. `PeriodicTable.__contains__`
  (`python/chargefw/_elements.py:58-62`) returns `False` on purpose. Add a three-line `__contains__`.

### A12. Smaller native API items

- **Status:** reported.
- `AtomParameterAccessor::operator[]` and `BondParameterAccessor::operator[]` delegate to checked
  `at()` (`src/parameters/models/parameter_view.cpp:13-15, 28-30`), unlike the unchecked `noexcept`
  indexing on parameter containers and classifications. Whether to align or document this difference
  remains a separate API decision; changing bounds semantics needs explicit approval.
- `include/chargefw/methods/method_prerequisites.h:32` and
  `include/chargefw/methods/method_applicability.h:36, 42-43`: `PrerequisiteIssue::kind`,
  `ExecutionIssue::kind`, and `ExecutionAssessment::mode` / `availability` have no default member
  initializer while their sibling members do, so `ExecutionAssessment a;` leaves the enums
  indeterminate. No misuse was found. No enumerator is a natural default. Resolution: leave; the
  library creates these values, and callers only read them.
- `src/calculation/assessment.cpp:251`: `to_string(MethodScope)` throws `std::logic_error` for an
  unknown value; the other `to_string` overloads in `calculation` throw `std::invalid_argument`.
- `include/chargefw/calculation/observer.h:14`: the public comment names the internal function
  `progress_for_indexed`.
- `include/chargefw/calculation/observer.h:98`: `CalculationCancelled` is public, but no public
  `calculate()` overload lets it escape (`src/calculation/calculation.cpp:219`). It could move to an
  internal header. See D11 for the related test-only seam.
- `src/core/atom.cpp:10` hard-codes `atomic_number > 100`, while
  `src/core/periodic_table.cpp:62-67` derives the bound from the table. Use
  `periodic_table().contains(atomic_number)`.

---

## D. Dead code

"No user" means that a repo-wide grep finds the name only at its declaration and definition. "Tests
only" means that the remaining users are under `tests/`. Re-run the grep before removing anything;
for public headers, removal is an API change and needs a CHANGELOG entry.

### D2. Four `CalculationInput` accessors — tests only

- **Status:** re-checked.
- **Where:** `include/chargefw/methods/calculation_input.h:61, 63, 67, 69`,
  `src/methods/calculation_input.cpp:53-60, 70-77`.
- **Check:** `grep -rln "has_geometry\|geometry_if_available\|has_parameters()\|parameters_if_available" include src apps python tests`
  lists the header, the `.cpp`, and `tests/cpp/methods/test_calculation_input.cpp`.
- **Fix:** remove the four and their test lines. Every method calls the throwing `geometry()` and
  `parameters()`.

### D3. `TopologyFeatures::bond_index_between` and `are_bonded` — tests only

- **Status:** re-checked.
- **Where:** `include/chargefw/features/topology_features.h:31-36`.
- **Check:** `grep -rln "bond_index_between\|are_bonded" include src apps python tests`.

### D4. `PdbReader::options()` and `MmcifReader::options()` — no user

- **Status:** reported.
- **Where:** `include/chargefw/adapters/gemmi/pdb_input.h:28`, `mmcif_input.h:31`.
  `PdbReader::options_` exists only for the accessor (`src/adapters/gemmi/pdb_input.cpp:224, 267`).
- **Check:** `grep -rn "\.options()" src apps python tests`.

### D5. Unused `ComplexityTerm` values and `to_string(ExecutionAvailability)`

- **Status:** re-checked for the two enumerators; reported for `to_string`.
- `ComplexityTerm::bonds` (`include/chargefw/methods/method_requirements.h:13`): referenced only in
  the notation switch at `src/methods/method.cpp:65`.
- `ComplexityTerm::bonds_squared` (`:17`): referenced only at `src/methods/method.cpp:73` and
  `src/methods/method_applicability.cpp:90`.
- `to_string(ExecutionAvailability)` (`include/chargefw/methods/method_applicability.h:26`): its only
  caller is `tests/cpp/methods/test_method_applicability.cpp:37`.
- **Check:** `grep -rn "ComplexityTerm::bonds\b\|ComplexityTerm::bonds_squared" include src apps python tests`.

### D6. Unreachable classification-count check

- **Status:** reported.
- **Where:** `src/methods/method_applicability.cpp:60-79` (`valid_classification_count`,
  `make_invalid_classification_count_issue`) and the `if` at `:265-271`.
- **Reason:** `:188-194` throws for methods with only common parameters, and the collection
  `check_parameter_prerequisites` returns either issues or exactly one classification per molecule
  (`src/methods/parameter_prerequisites.cpp:169-196`). `validate_selected_candidate` checks the count
  again at execution (`src/methods/applicable_method_execution.cpp:40-45`).

### D7. `breadth_first_distances(..., max_distance)` — parameter never used

- **Status:** reported.
- **Where:** `src/features/topology_helpers.h:9-11`, `src/features/topology_helpers.cpp:11-13, 26-28`.
- **Check:** `grep -rn "breadth_first_distances" include src apps python tests`; both callers
  (`topology_helpers.cpp:49, 93`) pass `std::nullopt`, and nothing outside the file calls it.
- **Fix:** drop the parameter and its branch; move the function into the anonymous namespace.

### D8. Line counter threaded through the MOL2 reader — never read

- **Status:** reported.
- **Where:** `src/adapters/native/mol2_input.cpp`: `std::size_t& line` is passed through `next_line`,
  `next_section`, `read_atoms`, `read_bonds`, `read_formal_charges`, and `parse_record`, but every
  MOL2 issue sets `.line = std::nullopt`.
- **Fix:** remove the parameter, or use it to report lines (which would be the better outcome).

### D9. `SourceIdentity.source_name` — no user

- **Status:** reported.
- **Where:** `python/chargefw/core.py:54-56`. An alias for `.source`.
- **Check:** `grep -rn "\.source_name" python tests docs`; every hit is on a `Molecule`.
  `docs/reference/python-api.md:68` lists only `source`, `record_index`, and `record_id`.

### D10. Public members used only by tests — maintainer's choice

- **Status:** reported.
- `PeriodicTable::contains(std::string_view)` (`include/chargefw/core/periodic_table.h:34`).
- `ParameterSet::publication()` and `notes()`
  (`include/chargefw/parameters/models/parameter_set.h:22-23`); `metadata()` exposes both.
- `ExecutionResult::cancelled()` (`include/chargefw/calculation/calculation.h:60`).
- `AtomicCharges::total()`, `ChargeSet::assignment(i)`, and `load_parameter_sets_json_directory` are
  also unused in production, but are reasonable conveniences for C++ consumers.
- `bundled_parameter_directory()` (`include/chargefw/parameters/io/parameter_data_paths.h:8`) is
  called only from `default_parameter_directories()` in the same `.cpp`.

### D11. `calculate(const CalculationRequest&)` and `CalculationResult` — tests only

- **Status:** reported.
- **Where:** `src/calculation/calculation_request.h:27-35`, `src/calculation/calculation.cpp:143-145`.
- **Reason:** production calls only the file-local `calculate_request` (`calculation.cpp:208`). The
  wrapper and the one-member `CalculationResult` are used by four test files.
- **Fix:** if kept as a test seam, return `charges::ChargeSet` directly and drop `CalculationResult`.
  The comment at `calculation_request.h:20` ("Zero delegates the worker count to the oneTBB runtime")
  disagrees with `thread_count()`, which uses `max_allowed_parallelism`.

### D13. Small dead code in `adapters`

- **Status:** reported.
- `src/adapters/gemmi/structure_import.h:25-29`: six default arguments of `make_record` are never
  used; both callers pass all eleven arguments. `if (name.empty()) name = structure.name;` at
  `structure_import.cpp:350` cannot change anything given what both callers pass.

---

## C. Duplicated code

Each item is meant to become one helper. Larger merges that the review considered and rejected are
listed at the end of the section.

### C1. Charge-constrained dense solve and its cost estimate

- **Status:** reported.
- **Where (solve):** `src/methods/builtin/eem.cpp:82-85, 99-102`, `sfkeem.cpp:78-85`,
  `qeq.cpp:112-119`, `eqeq.cpp:92-99`, `tsef.cpp:83-90`, `smpqeq.cpp:94-99`, `abeem.cpp:192-197`.

  ```cpp
  A.row(n).setOnes();
  A.col(n).setOnes();
  A(n, n) = 0.0;
  b(n) = input.target_charge();
  const Eigen::VectorXd q = solve_in_place(A, b).head(n);
  return charges::AtomicCharges{std::vector<double>{q.data(), q.data() + q.size()}};
  ```
- **Where (estimate):** byte-identical in `eem.cpp:105-112`, `sfkeem.cpp:88-95`, `qeq.cpp:122-129`,
  `eqeq.cpp:106-113`; variants in `smpqeq.cpp:111-118` and `eqeqc.cpp:83-91`.
- **Fix:** two inline helpers in `src/methods/builtin/dense_solve.h`:
  `solve_charge_constrained(matrix, rhs, target_charge)` and `charge_constrained_estimate(n)`. Leave
  `src/parameterization/eem_derivatives.cpp`, which keeps its own LU on purpose.
- **Validation:** numerical code; run the method tests in a release profile as well.

### C2. Copying one atom parameter into a vector

- **Status:** reported.
- **Where:** `src/methods/builtin/sfkeem.cpp:48-52`, `smpqeq.cpp:55-59`, `qeq.cpp:75-79`,
  `eqeqc.cpp:48-52`, `sqe.cpp:71-75`, `sqeqp.cpp:29-33`,
  `src/parameterization/sqe_derivatives.cpp:79-83`.
- **Reason it exists:** `AtomParameterAccessor::operator[]` does a classification lookup and a by-name
  parameter lookup on every call (`src/parameters/models/parameter_view.cpp:13-20`).
- **Fix:** one helper `atom_values(view, name, atom_count) -> std::vector<double>` in a small internal
  header under `src/methods`.

### C4. Parameter-set JSON reader

- **Status:** reported.
- `src/parameters/io/parameter_set_io.cpp:332-372` (`parse_atom_parameters`) and `:374-414`
  (`parse_bond_parameters`) differ in the section name, the key parser, and the entry type. Fix: one
  function template taking the key parser.
- `:486-491` and `:693-698` repeat the directory loading loop. Fix: extract
  `append_parameter_sets_from_directory(sets, ids, directory)`. The `_if_present` variant checks
  `is_directory` at `:481-484`, and `json_files_in_directory` checks it again at `:441`.
- `src/parameters/models/parameter_set.cpp:11-16` (`same_names`) and the lambda at
  `parameter_set_io.cpp:592-594` are the same comparison. Fix: move `same_names` into
  `src/parameters/parameter_helpers.h`.

### C6. PEOE, MPEOE, and GDAC helpers

- **Status:** reported.
- `uses_formal_initial_charges()`: `src/methods/builtin/peoe.cpp:26-34`, `mpeoe.cpp:26-34`.
- `is_hydrogen()`: `peoe.cpp:15-17`, `mpeoe.cpp:16-18`.
- `atom_denominator()`: `mpeoe.cpp:20-24`, `gdac.cpp:19-23`.
- Formal-charge initialization: `peoe.cpp:58-63`, `mpeoe.cpp:59-64`, and a third time in
  `sqe_core::formal_initial_charges` (`sqeq0.cpp:13-23`).
- Within `peoe.cpp` and `mpeoe.cpp`, the prerequisite reads the `initial_charges` option through the
  helper, while `calculate()` reads it with `get<std::string>("initial_charges") == "formal"`.
- **Related behaviour difference:** GDAC throws on a zero denominator (`gdac.cpp:75`); PEOE and MPEOE
  divide unchecked and fail later with the generic "atomic charges must be finite".
- **Fix:** move the shared helpers into `prerequisite_helpers.h` or a small `peoe_common.h`, and use
  the option helper in `calculate()` too.

### C7. Derivative code repeats checks of the methods

- **Status:** reported.
- `require_finite` exists three times: `src/parameterization/eem_derivatives.cpp:112-118`,
  `sqe_derivatives.cpp:223-229`, `normal_equations.cpp:120-124`.
- `sqe_derivatives.cpp:254-258` re-implements the validation at the top of `sqe_core::calculate`
  (`src/methods/builtin/sqe.cpp:183, 193-195`). The two are documented to fail identically and are
  kept in step by hand.
- `SqeqpLinearization::jacobian` (`sqe_derivatives.cpp:303-316`) copies the predicate of
  `check_columns` (`src/parameterization/parameter_columns.h:22-28`).
- **Fix:** one `require_finite(std::span<const double>, std::string_view what)` in
  `parameter_columns.h`; expose the SQE preamble as `sqe_core::validate_input(input, initial)` in
  `sqe_system.h`; call `check_columns` in the SQE+qp jacobian.

### C8. Prefixing an issue with its molecule

- **Status:** reported.
- **Where:** `src/methods/method_prerequisites.cpp:12-21` (`with_molecule_prefix`) and inline in
  `src/methods/parameter_prerequisites.cpp:176-183`.
- **Fix:** move the function to `src/methods/parameter_issue_messages.h` and use it in both.

### C9. Optional parameter-set ID built by hand

- **Status:** plausible. The helper tests `uses_parameters()`, the copies test
  `parameter_set == nullptr`; check that they are equivalent (`ApplicableMethod::uses_parameters()`
  is defined as `parameter_set != nullptr`).
- **Where:** `src/calculation/assessment.cpp:373-375`, `src/calculation/calculation.cpp:191-193`
  repeat `methods::detail::parameter_set_id_for` (`src/methods/applicable_method_execution.cpp:63`).

### C10. Registry lookup with an "is not registered" error

- **Status:** reported.
- **Where:** `src/calculation/assessment.cpp:83-86, 103-106, 152-155`;
  `apps/chargefw/reports.cpp:225-228, 278-280`; `apps/chargefw/input.cpp:69-72` (other wording).
- **Fix:** one `require_method(id) -> const Method&` next to the registry, or a file-local helper in
  `assessment.cpp`.

### C11. Duplicates in the CLI

- **Status:** reported.
- `apps/chargefw/input.cpp:133-136, 141-144, 149-152, 157-160`: the same "Structural input options
  are only supported for PDB and mmCIF input" check four times. Hoist it before the dispatch.
- `apps/chargefw/output.cpp:202-206` and `:222-227`: the same metrics, `write_json`, and
  `report_diagnostics` sequence in the no-charges and success branches. The peak-memory sample is
  also taken twice (`:188`, `:221`).

### C12. Output writers

- **Status:** reported.
- `safe_line`: identical in `src/adapters/native/mol2_output.cpp:24-28` and `pqr_output.cpp:25-29`.
- The "Generated by ChargeFW …; method=…; parameter_set=…" line: `mol2_output.cpp:125-133`,
  `pqr_output.cpp:201-209`.
- Record name fallback (record ID, molecule name, `molecule_N`): `mol2_output.cpp:33-37`,
  `src/adapters/gemmi/mmcif_output.cpp:94-98`.
- **Fix:** `generator_comment(charge_set, version)` and `record_display_name(record, index)` in an
  internal adapters header.

### C13. Input readers

- **Status:** reported.
- Reading a line and stripping `\r`: `src/adapters/native/common_input.cpp:14-28`,
  `mol2_input.cpp:45-55, 329-332`, `sdf_input.cpp:33-36`, `src/adapters/gemmi/pdb_input.cpp:110-113`.
  The MOL2 and SDF readers bypass `common_input` only because `read_line` throws at the end of input.
  Fix: add `common::try_read_line(input) -> std::optional<std::string>` and build `read_line` on it.
- "Each atom is in exactly one component": `src/adapters/fixed_ions.cpp:118-148` and
  `src/adapters/gemmi/mmcif_output.cpp:346-368`. Fix: one `detail::component_of_atom(record)` in
  `src/adapters/component.h`.
- Non-finite coordinate test: `src/features/conformer_features.cpp:41-44` and
  `src/features/spatial_index.h:21-28`. Fix: call `geometry.first_nonfinite_atom_index()` in the
  point-cloud constructor.

### C14. Python package

- **Status:** reported.
- "Reject bool, `operator.index`, non-negative" in six places:
  `python/chargefw/_calculation_options.py:21-32`, `charges.py:35-44`, `core.py:40-47`,
  `_calculation_values.py:480-505`, `io/rdkit.py:169-185`, `parameterization/_evaluator.py:265-270`.
  They already disagree: `rdkit.attach_charges` raises `ValueError` for a negative conformer where
  `CalculationResult.assignment` raises `IndexError`, and two of them reject only `bool`, not
  `np.bool`. Fix: one helper in `_arrays.py`.
- Wrapping a one-set native catalog as a `ParameterSet`: the bundled-ID check
  (`calculation.py:190-194`, `parameterization/_template.py:31-32`, `_fit.py:76-77`), the
  descriptor unwrap (`calculation.py:189, 195`, `_template.py:38-39`, `_fit.py:82-83`), and
  `parameter_set._native or _default_parameter_catalog()` (`calculation.py:205`, `_fit.py:78`,
  `_evaluator.py:327`). Fix: two private helpers in `calculation.py`.
- `io/rdkit.py:137-150` (`_assignment`) picks an assignment by other rules than
  `CalculationResult.assignment` (`_calculation_values.py:477-511`): it falls back silently where the
  latter raises. Plausible; check the RDKit tests before changing error types.
- `parameterization/_diagnostics.py:97-105` repeats `_compression.py:51-63` (scale to unit diagonal,
  `eigh`, undetermined check). The first clips a negative diagonal, the second raises. Low value.

### C15. Python extension

- **Status:** reported.
- Observer setup: the same seven lines in `python/src/calculation.cpp:337-342` and `:465-470`.
- `method_option_value` (`python/src/calculation.cpp:145-147`) is identical to `option_value`
  (`python/src/methods.cpp:17-19`).

### Considered and not recommended

- Turning the `Atom`/`Bond` class pairs (`AtomParameters`/`BondParameters`, the two classifications,
  the two accessors) into class templates. It removes about 200 lines but changes the public types.
- Merging the PEOE, MPEOE, and GDAC iteration loops. They differ in damping and denominator, and
  each must stay faithful to its cited method.
- Moving the "selection and bond options only for PDB and mmCIF" rule and the `hybrid` default from
  the CLI (`apps/chargefw/input.cpp:133-160`) and Python (`python/chargefw/io/__init__.py:143-159,
  195`) into the native adapters. It is real policy duplication, but fixing it is a native API change.

---

## X. Checked and found clean, or rejected

Listed so that these are not raised again without new evidence.

- **X1. Fixed ion at the position of a calculated atom.** One reader suspected an undiagnosed
  division by zero in `src/methods/builtin/eem.cpp:31-32` and `sqe.cpp:101-102`. Rejected:
  `src/calculation/fixed_charge_partition.cpp:110-129` rejects a source that coincides with an active
  atom before any method runs.
- **X2. NaN charges in mmCIF output.** `src/adapters/gemmi/mmcif_output.cpp:670` tests
  `charge < -10.0 || charge > 10.0`, which NaN passes. Not reachable through the library, because
  `charges::AtomicCharges` rejects non-finite values at construction. Writing
  `!(charge >= -10.0 && charge <= 10.0)` would make the check self-sufficient.
- **X3. Eigen aliasing in DENR.** `q = lu.solve(q - tmp);` in `src/methods/builtin/denr.cpp:63` was
  probed with a standalone program against the pinned Eigen, with a matrix that forces a non-identity
  permutation; results matched the `.eval()` form at `-O0` and `-O2`.
- **X4. Parallel execution.** `src/calculation/parallel_for.h`: exceptions propagate through
  `tbb::parallel_for`, which joins all tasks before rethrowing, so by-reference captures of locals
  are safe. Fragment progress is serialized by a mutex; observer exceptions are swallowed. Full
  execution parallelizes targets and reduced execution parallelizes fragments; they do not nest.
- **X5. Global mutable state.** None in the library. The registry, derivative tables, method
  metadata, and option schemas are function-local `static const`. The only mutable static is the
  CLI's interrupt flag.
- **X6. Borrowed-pointer types.** `TopologyFeatures`, `ConformerFeatures`, `PreparedMolecule`,
  `ParameterView`, `SpatialFragmentBuilder`, and `CalculationInput` delete their rvalue overloads
  consistently. `AssessmentResult` keeps molecules, partition, and prepared molecules on the heap and
  deletes move assignment, so plan pointers survive move construction.
- **X7. Python bindings.** Every array returned by the extension owns a copy; no view of native memory
  is returned. Plans and evaluators share the assessment by `shared_ptr`. Observer callbacks take the
  GIL, and the observer is destroyed after the GIL is re-acquired. Input arrays are declared
  C-contiguous with fixed dtype and rank. The `.pyi` stubs match the registrations.
- **X8. Parsers.** Bond endpoints are range-checked in every reader or by the `core::Molecule`
  constructor. Integer parsing uses `from_chars` with a full-consumption check. Output uses
  `std::format` without a locale. mmCIF free text goes through `gemmi::cif::quote`.
- **X9. Not examined further.** Duplicate PDB serials make a CONECT record resolve to the first atom
  with that serial (`src/adapters/gemmi/selection.cpp:64`, `bonds.cpp:357-368`); what Gemmi does with
  overflowed serials was not checked. EQeq+C accepts elements with a zero covalent radius (Bk, Cf,
  Es, Fm), which ABEEM rejects (`src/methods/builtin/abeem.cpp:62-64`); it matters only if a
  parameter set supplies those elements. The library path from the dynamic loader may be relative
  (`src/platform/linux/library_location.cpp:24`); it fails only if the process changes directory
  before loading parameters. A right-trimmed V2000 atom line of 32 characters is rejected
  (`src/adapters/native/mol_input.cpp:117`), although the comment at `:46` says writers may end the
  line after the symbol.

  If the PDB serial or the trimmed V2000 line is confirmed, a sentence in
  `docs/reference/formats.md` is the likely resolution, because both come from unusual files: say
  that `CONECT` records need unique serials, and that V2000 atom lines need the full fixed width.
  The other two need neither a fix nor a notice unless someone reports them.
