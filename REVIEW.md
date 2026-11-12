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
| X9 (PDB serials, trimmed V2000 lines) | Document if confirmed | `docs/reference/formats.md` |

Follow [AGENTS.md](AGENTS.md) when acting on a finding. In particular, the findings marked
*design-sensitive* change public API or layering and need a sketch and the maintainer's review before
implementation. Remove a finding from this file when it is fixed or rejected; do not keep history
here.

The greps below are meant to be run from the repository root. "Repo-wide" means
`include src apps python tests`.

## Summary

No memory-safety bug and no data race was found in the library.

| Section | Content | Findings |
| --- | --- | --- |
| B | Behaviour bugs | No remaining findings |
| A | Public API design | No remaining findings |
| D | Dead code | D4–D13 |
| C | Duplicated code | C1–C15 |
| T | Tests that protect nothing | T1–T25 |
| X | Checked and found clean, or rejected | X1–X9 |

---

## D. Dead code

"No user" means that a repo-wide grep finds the name only at its declaration and definition. "Tests
only" means that the remaining users are under `tests/`. Re-run the grep before removing anything;
for public headers, removal is an API change and needs a CHANGELOG entry.

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
