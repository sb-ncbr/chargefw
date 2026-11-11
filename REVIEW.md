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
- **Line numbers refer to commit `cff64e0`**, except in section T, which states its own. If they no
  longer match, search for the quoted code.
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
| B | Behaviour bugs | B2–B8 |
| L | Lifetime and misuse traps | L1–L2 |
| A | Public API design | A4–A12 |
| D | Dead code | D4–D13 |
| C | Duplicated code | C1–C15 |
| T | Tests that protect nothing | T1–T25 |
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

### B6. Unsupported V3000 blocks need a representation decision

- **Status:** re-checked for rejection of blocks other than `ATOM`, `BOND`, and `COLLECTION`.
- **Where:** `src/adapters/native/mol_input.cpp` (`parse_v3000`).
- **Problem:** valid records containing `SGROUP` or `3D` blocks are rejected. Some S-groups carry
  query or bond semantics, so skipping every unknown block can change the interpreted structure.
- **Resolution (design-sensitive):** determine which additional blocks can be ignored safely and
  which require interpretation or explicit rejection. Do not replace the current rejection with
  blanket block skipping. The supported subset is documented in `docs/reference/formats.md`.

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
  keeps its molecule pointer while its vectors are emptied, so `degree`, `neighbor_indices`, and
  `incident_bond_indices` pass `validate_atom_index` (which checks the
  molecule's atom count) and then index an empty vector.
- **Fix:** validate against `neighbor_indices_.size()`. The alternative is to copy `topology_` in the
  move operations, as is done for the cache.
- **Resolution:** document. No workflow uses a moved-from prepared molecule. Correct the comment
  above the move constructor in `prepared_molecule.h`: only `geometry_defects()` stays valid on the
  moved-from object; its topology must not be queried. This is not user documentation.

---

## A. Public API design

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

### A7. `JsonWriter::write` does not check stream failure

- **Status:** reported.
- **Where:** `include/chargefw/adapters/native/json_output.h:16`.
- **Problem:** It is the only writer that does not throw when the stream fails; the others end with
  `if (!*output_) throw`. The CLI checks the stream itself, so only direct C++ users are affected.
- **Fix:** add the same stream check.
- **Resolution:** either add the check or say in the comment at the declaration that the caller checks
  the stream after `write()`; the example in `docs/interfaces/cpp.md:162` should then show that check.

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
- `include/chargefw/calculation/observer.h:98`: `CalculationCancelled` is public, but no public
  `calculate()` overload lets it escape (`src/calculation/calculation.cpp:219`). It could move to an
  internal header. The internal `calculate(CalculationRequest)` overload lets it propagate.

---

## D. Dead code

"No user" means that a repo-wide grep finds the name only at its declaration and definition. "Tests
only" means that the remaining users are under `tests/`. Re-run the grep before removing anything;
for public headers, removal is an API change and needs a CHANGELOG entry.

### D5. Unused `ComplexityTerm` values and `to_string(ExecutionAvailability)`

- **Status:** re-checked for the two enumerators; reported for `to_string`.
- `ComplexityTerm::bonds` (`include/chargefw/methods/method_requirements.h:13`): referenced only in
  the notation switch at `src/methods/method.cpp:65`.
- `ComplexityTerm::bonds_squared` (`:17`): referenced only at `src/methods/method.cpp:73` and
  `src/methods/method_applicability.cpp:90`.
- `to_string(ExecutionAvailability)` (`include/chargefw/methods/method_applicability.h:26`): no callers.
- **Check:** `grep -rn "ComplexityTerm::bonds\b\|ComplexityTerm::bonds_squared" include src apps python tests`.

---

## C. Duplicated code

Each item is meant to become one helper. Larger merges that the review considered and rejected are
listed at the end of the section.

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

### C10. Registry lookup with an "is not registered" error

- **Status:** re-checked.
- **Where:** `apps/chargefw/reports.cpp:225-228, 278-280` and
  `python/src/parameterization.cpp:241-244`.
- **Problem:** these repeat the same registry lookup and `std::invalid_argument` for an unknown ID.
- **Option (design-sensitive):** add `MethodRegistry::at(id) -> const Method&` and use it here and
  in assessment. This public API addition needs a separate decision.
- Keep `apps/chargefw/input.cpp:69-72` separate: its `std::runtime_error` and message describe a CLI
  option-parsing error.

### C12. Output writers

- **Status:** reported.
- `safe_line`: identical in `src/adapters/native/mol2_output.cpp:24-28` and `pqr_output.cpp:25-29`.
- The "Generated by ChargeFW …; method=…; parameter_set=…" line: `mol2_output.cpp:125-133`,
  `pqr_output.cpp:201-209`.
- Record name fallback (record ID, molecule name, `molecule_N`): `mol2_output.cpp:33-37`,
  `src/adapters/gemmi/mmcif_output.cpp:94-98`.
- **Fix:** `generator_comment(charge_set, version)` and `record_display_name(record, index)` in an
  internal adapters header.

### C14. Python package

- **Status:** remaining behavior differences need a separate decision.
- `rdkit.attach_charges` raises `ValueError` for a negative conformer while
  `CalculationResult.assignment` raises `IndexError`. Range checks remain local to preserve these
  exception types.
- `io/rdkit.py:137-150` (`_assignment`) picks an assignment by other rules than
  `CalculationResult.assignment` (`_calculation_values.py:477-511`): it falls back silently where the
  latter raises. Plausible; check the RDKit tests before changing error types.
- `parameterization/_diagnostics.py:97-105` repeats `_compression.py:51-63` (scale to unit diagonal,
  `eigh`, undetermined check). The first clips a negative diagonal, the second raises. Low value.

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
