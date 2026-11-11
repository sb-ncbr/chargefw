# Refactor: execution plans own their assessment data

Handoff notes for continuing this work. Base commit: `c41320c`. The change is in the working tree,
not staged or committed. Read `AGENTS.md` first.

## Goal

An `ExecutionPlan` used to hold raw pointers into its `AssessmentResult`, so a plan was only valid while
its assessment lived, and every plan-taking function needed both objects plus a runtime check that they
matched. Now a plan shares ownership of the assessed data and is complete on its own.

## What changed

### Ownership

- New private `detail::AssessmentState` (`src/calculation/assessment_state.h`): molecules, parameter
  sets, fixed-ion partition, prepared molecules, applicability result, rejections, input warnings, and
  assessment time. Created once on the heap, never moved (members refer to each other), filled during
  `AssessmentResult::assess_owned()` and then only read.
- `AssessmentResult` holds `std::shared_ptr<const detail::AssessmentState>` plus `std::vector<ExecutionPlan>`.
  The plans live outside the state, so there is no ownership cycle.
- `ExecutionPlan` holds the same `shared_ptr` (replacing the empty `PlanIdentity` tag, now deleted),
  a pointer to its candidate inside the state, its policy, and its warnings. It is copyable; a copy
  keeps the whole assessment state alive.

### Public C++ API

| Before | After |
| --- | --- |
| `calculate(assessment, plan, threads, observer)` | `calculate(plan, threads, observer)` |
| `calculate(assessment, plan, candidate, ...)` | `calculate(plan, candidate, ...)` |
| `estimate(assessment, plan, threads)` | `estimate(plan, threads)` |
| `parameterization::normal_equations(assessment, plan, ...)` | `normal_equations(plan, ...)` |
| `entry_counts`, `atom_entries`, `single_bond_entries(assessment, plan)` | `(plan)` |
| `ExecutionPlan` move-only | copyable |
| `plans()`, `default_plan()` callable on a temporary | `const&&` overloads deleted |

- `calculate(assessment)` (default plan) is unchanged.
- `AssessmentResult` stays move-constructible and not move-assignable, so pointers to its plans stay
  valid while it exists.
- The "plan belongs to a different assessment" `std::invalid_argument` no longer exists, and the two
  tests for it were removed.
- No deprecated forwarding overloads were added. Outside C++ callers must drop the first argument.

### Internal

- `ExecutionPlan::execute()` and `ExecutionPlan::with_values()` replace
  `AssessmentResult::execute()`, `check_plan()`, and `replace_values()` (`src/calculation/calculation.cpp`).
- `detail::AssessmentAccess` (`src/calculation/assessment_access.h`) now reads through a plan:
  `source_molecules(plan)`, `molecules(plan)` (prepared), `partition(plan)`, `classifications(plan)`,
  `with_values(plan, candidate)`.
- Callers updated: `apps/chargefw/{main,reports}.cpp`, `python/src/{calculation,parameterization}.cpp`
  (binding internals only; the Python API is unchanged), and the C++ tests.

### Documentation

- `docs/interfaces/cpp.md`: new call forms; `plans()`/`default_plan()` lifetime; a copied plan keeps the
  molecules and parameter sets in memory; an assessment can be moved but not assigned.
- `DEVELOPMENT.md`: the parameterization note on reading assessment data through a plan.

## Validation done

- `gcc-debug`: builds without warnings; all 39 native tests pass.
- Python tests in a temporary uv venv: 126 pass. `test_gemmi.py`, `test_output.py`, and
  `test_examples.py` were not run: `gemmi==0.7.4` has no wheel for the system Python 3.15.
- `clang-asan`: `test_calculation`, `test_fixed_ion_planning`, `test_normal_equations`, and
  `test_parameter_entries` pass, including a test that calculates with a plan copied out of an
  assessment that was already destroyed.

Not run: `clang-debug`, `clang-ubsan`, `clang-tidy`, mypy with gemmi, the full ASan suite.

## Machine constraints

- Never build with Ninja's default parallelism: it triggered the OOM killer. Use `-j4` or lower for
  debug and release profiles, and `--parallel 1` for `clang-asan` and `clang-ubsan`. Never run two
  builds at once.
- Python test dependencies missing from the system Python go into a temporary uv venv using
  `/usr/bin/python3`; run pytest from `build/<preset>/tests/python` with
  `PYTHONPATH=build/<preset>/python`.

## Open follow-ups

1. **Python API.** `chargefw.calculate(molecules, plan)` takes the molecules only to check they are
   the same objects (`python/chargefw/calculation.py`, `Plan._matches`). Mirror the C++ change with
   `plan.calculate()` or `chargefw.calculate(plan)`. The binding's `NativePlan` keeps
   `shared_ptr<NativeAssessmentState>` plus an index; that state is still needed for the imported
   records and the request, so holding a copied `ExecutionPlan` there gains little.
2. **Reusing a plan's choices on other molecules.** Proposed, not designed: `Plan.settings`, a fully
   explicit `RequestedCalculation` (method, parameter set, options, parameter matching, execution mode
   and radius, fixed ions), accepted by `calculate(molecules, settings)` and `assess()`. It reassesses
   the new molecules with exactly those choices and fails instead of falling back. The name is open:
   "requested" describes what was asked, not reusable settings. Design-sensitive: sketch it and get
   user review first.
3. **Deprecation.** Decide whether outside C++ users of 3.0.0 need `[[deprecated]]` overloads with the
   old `(assessment, plan, ...)` signatures for one release.
4. **Remaining qualification** listed under "Not run" above, before release.
