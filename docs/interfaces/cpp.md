# Native C++ library

ChargeFW exposes a C++23 shared library through the `chargefw::core` CMake target. The public API owns
molecules, parameter sets, assessments, plans, results, and charge arrays; ordinary application code does
not need to manage prepared-feature or classification lifetimes.

## Namespaces

The public namespaces group the API by responsibility:

| Namespace | Responsibility |
| --- | --- |
| `chargefw::core` | Owned atoms, bonds, molecules, conformers, and periodic-table data |
| `chargefw::features` | Prepared topology, geometry, and spatial fragments |
| `chargefw::parameters` | Parameter data, loading, and atom/bond classification |
| `chargefw::methods` | Method metadata, options, requirements, applicability, and algorithms |
| `chargefw::calculation` | Assessment, deterministic planning, execution policy, progress, and calculation |
| `chargefw::charges` | Source-indexed charge assignments and charge sets |
| `chargefw::adapters` | Native and Gemmi-backed molecular input and output |

## Build and link

ChargeFW requires CMake 3.28 or newer, Ninja, a GCC or Clang toolchain with C++23 support, and installed
development packages for nlohmann/json 3.11.3 or newer, Eigen 3.4 or 5, nanoflann 1.7 or newer,
oneTBB 2022.1 or newer, Gemmi 0.7.4 or newer, and, for the command-line application, CLI11 2.4 or newer.
To build without them, set `CHARGEFW_FETCH_DEPENDENCIES=ON`; ChargeFW then downloads pinned,
checksum-verified sources of all of them.

```bash
cmake -S . -B build/native-release -G Ninja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INTERPROCEDURAL_OPTIMIZATION=ON \
    -DCHARGEFW_BUILD_CLI=OFF \
    -DCHARGEFW_BUILD_PYTHON=OFF \
    -DCHARGEFW_BUILD_TESTS=OFF \
    -DCHARGEFW_ENABLE_NATIVE_OPTIMIZATIONS=OFF
cmake --build build/native-release --parallel
cmake --install build/native-release --prefix "$PWD/_install" --strip
```

Set `CHARGEFW_ENABLE_NATIVE_OPTIMIZATIONS=ON` only for an installation that will remain on the build
machine. IPO does not make the installed binaries CPU-specific.

Use the installed CMake package from another project:

```cmake
find_package(chargefw CONFIG REQUIRED)
target_link_libraries(my_application PRIVATE chargefw::core)
```

The installed package requires the exact nlohmann/json and Gemmi versions that ChargeFW was built
against, because their types appear in the public headers. With installed dependency packages,
downstream CMake projects need those packages discoverable and the oneTBB runtime library available.
With `CHARGEFW_FETCH_DEPENDENCIES=ON`, the installation also includes the nlohmann/json and Gemmi
development packages and the oneTBB runtime.

The installed prefix can be relocated as a unit. `load_default_parameter_sets()` resolves bundled JSON
relative to the loaded ChargeFW library.

## Basic calculation

```cpp
#include <chargefw/calculation/calculation.h>
#include <chargefw/core/atom.h>
#include <chargefw/core/bond.h>
#include <chargefw/core/conformer.h>
#include <chargefw/core/molecule.h>
#include <chargefw/core/molecule_collection.h>
#include <chargefw/parameters/io/parameter_set_io.h>

#include <iostream>
#include <utility>
#include <vector>

int main() {
    using namespace chargefw;

    auto atoms = std::vector<core::Atom>{core::Atom{8}, core::Atom{1}, core::Atom{1}};
    auto bonds = std::vector<core::Bond>{core::Bond{0, 1}, core::Bond{0, 2}};
    auto conformers = std::vector<core::Conformer>{core::Conformer{{
        {0.0, 0.0, 0.0},
        {0.96, 0.0, 0.0},
        {-0.24, 0.93, 0.0},
    }}};
    auto molecule = core::Molecule{
        std::move(atoms), std::move(bonds), std::move(conformers), "water"};

    auto request = calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{
            std::vector<core::Molecule>{std::move(molecule)}},
        .parameter_sets = parameters::load_default_parameter_sets(),
        .method_id = "eem",
        .execution_selection = calculation::ExecutionSelection{
            calculation::ExecutionSelectionKind::full},
    };

    auto assessment = calculation::assess(std::move(request));
    auto result = calculation::calculate(assessment);
    if (!result.calculated()) {
        std::cerr << result.failure_message.value_or("calculation failed") << '\n';
        return 1;
    }

    for (double value : result.charges->assignment(0).charges.values()) {
        std::cout << value << '\n';
    }
}
```

Passing an lvalue `AssessmentRequest` copies its molecule and parameter inputs. Passing an rvalue
transfers them into the assessment. An `AssessmentResult` owns the prepared molecules, classifications,
and plans needed for execution.

## Molecule model

`chargefw::core::Molecule` owns:

- a source-ordered vector of `Atom` values with atomic number, formal charge, and optional source name;
- indexed single, double, or triple `Bond` values;
- zero or more `Conformer` values, each containing one position per atom; and
- an optional molecule name.

`MoleculeCollection` owns source-ordered molecules and an optional collection name. Molecules accept the
atomic numbers 1–100 represented by the bundled periodic table; individual methods and parameter sets may
support a smaller subset. Constructors validate atomic numbers, bond endpoints and duplicates, conformer
cardinality. Non-finite coordinates are retained, but methods that require geometry are inapplicable when
coordinates are missing or non-finite.

Lower-level `methods::CalculationInput` accepts a trailing `std::span<const methods::FixedPointSource>`;
each source contains a Cartesian `core::Position` and charge. The input borrows the span, so keep its
elements alive for the input's use. These positional sources are distinct from the original-collection
indices in `calculation::FixedIons`; direct method callers provide positions themselves, while
the assessment facade validates indexed selectors and constructs sources for methods that support
fixed point sources.
`MethodRequirements::supports_fixed_point_sources` reports numerical support. Direct callers of
`Method::calculate()` must satisfy the method's requirements, including this capability; methods that do
not support fixed point sources require an empty span. The flag describes a caller precondition;
it does not intercept virtual calls. The source coupling is defined by the selected method; the facade's
effective provenance already records the method, parameter set, and method options.

EEM, SQE, SQE+q0, and SQE+qp support fixed point sources in full, cutoff, and cover facade execution.
Their facade `CalculationInput::target_charge()` is set to the prepared active formal-charge sum, and
direct method results cover only active atoms; direct callers supply method inputs and handle any further
assembly. Reduced fragments are built from active atoms and receive the complete fixed-source field,
including sources outside the fragment radius. Conservation corrections apply to active charges; the
facade restores fixed values once in original atom order. SQE requires neutral active connected
components, SQE+q0 retains their prepared formal-charge seeds, and SQE+qp normalizes its fitted reference
charges globally to the active target before bond-transfer equilibration.
The [Fixed ions](../guide/fixed-ions.md) describes active-charge provenance.
Direct callers must supply finite source values and finite geometry for the active molecule, and ensure
that source positions do not coincide with active atom positions. The
[point-source terms](../guide/fixed-ions.md#point-source-terms) describes the point-source field term.

## Assessment

`calculation::AssessmentRequest` contains:

- the owned `MoleculeCollection` and parameter sets;
- an optional method ID and its parameter-set ID;
- method-scoped `methods::MethodOptions` overrides;
- strict or permissive parameter classification options;
- an `ExecutionSelection`;
- a `ResourcePolicy` containing automatic execution thresholds; and
- an optional `FixedIons`, whose sources pair a fixed charge with zero-based molecule and
  atom indices in the supplied collection.

For ordinary calculations, leave `fixed_ions` absent or provide an empty `sources` vector.
Native callers prescribe indexed charges directly, independently of the CLI/Python ion catalog.
Nonempty sources must use distinct in-range molecule/atom indices, finite charges, atoms with no graph
bonds, and leave at least one active atom in each affected molecule. Each affected molecule needs
conformer coordinates that are finite throughout; source coordinates must not exactly coincide with an
active atom in the same conformer. The active charge total is the sum of formal charges on the unselected
atoms in the prepared molecule; selected atoms' imported formal charges do not alter that active total.
The modeled molecular total is the active total plus the prescribed source charges.
A valid nonempty selection filters methods
by their fixed-ion capability and assesses execution modes. EEM, SQE, SQE+q0, and SQE+qp can produce
full, cutoff, and cover plans. Automatic selection may choose an un-warned reduced plan when full carries
a resource warning; explicit full selection may execute with its warning. Parameter classification and
execution use the active molecules. Successful facade results restore fixed charges in original atom
order and expose fixed atom indices and assigned charges through `ExecutionResult::effective`. The
[Fixed ions](../guide/fixed-ions.md) describes this boundary.
Invalid source selections still throw `std::invalid_argument`.

`calculation::assess()` returns an `AssessmentResult` with read-only access to its original owned source
`molecules()`, priority-ordered `plans()`, structured `rejections()`, `default_plan()`, and applicability
timing. A plan exposes its method, parameter set, validated method options, concrete `ExecutionPolicy`,
and warnings. Plans are tied to the assessment that created them. `calculation::estimate(assessment,
plan, max_threads)` returns the plan's approximate peak working memory and operation count as a
`methods::ResourceEstimate`; the [resource estimates](../guide/large-structures.md#threads-and-estimates) describes the
model.

Use `methods::method_registry()` to inspect built-in method metadata, human-readable implementation
notes, and option schemas. Use `methods::complexity_notation()` to format declared time and memory
complexity, and `Method::estimate()` for the approximate memory and work of one target of a given
`TargetSize` with validated options. The notation, variables, and estimate accounting are defined in the
[Methods and parameter sets](../guide/methods.md). Use
`parameters::load_default_parameter_sets()` for the installed bundled catalog. Parameter-set IDs in one
assessment request must be unique. The [parameter-set JSON reference](../reference/parameter-format.md) documents the bundled
data format and classifier semantics.

## Execution policy

`ExecutionSelectionKind` accepts automatic, full, cutoff, or cover selection. Effective plans always use
the concrete `ExecutionMode` values full, cutoff, or cover. `std::nullopt` disables an automatic resource
threshold. The [Large structures](../guide/large-structures.md) defines mode semantics, defaults,
supported methods, conservation behavior, and approximation limits.

Pass the thread limit to `calculate()` for each execution. A count of zero delegates scheduling to
oneTBB, and the same prepared assessment may be executed repeatedly with different limits:

```cpp
auto assessment = calculation::assess(std::move(request));
auto serial = calculation::calculate(assessment, 1);
auto parallel = calculation::calculate(assessment, 4);
```

Method implementations declare reduced-execution support through
`ResourceRequirements::reduced_charge_policy` and `ReducedChargePolicy`. Each plan separately exposes its
selected execution mode and radius through its concrete `ExecutionPolicy`.

## Results, errors, and mapping

`calculation::ExecutionResult` reports an `ExecutionStatus`, optional `charges::ChargeSet`, rejections,
effective calculation provenance, optional failure text, and applicability/computation timings.

`adapters::ChargeCalculationResult` is the immutable application-facing boundary. It owns the ordered
`ImportedMoleculeRecord` values used to build the calculation request, requested provenance, and the
`ExecutionResult`. JSON and generated molecular writers can therefore consume the result without a
separately supplied molecule collection, identity list, diagnostics list, or import context.

Result construction validates import-mapping dimensions, successful assignment coverage, canonical
molecule/conformer order, uniform scope, atom dimensions, target bounds, and agreement with effective
method and parameter-set provenance. Optional `EffectiveCalculation::fixed_ions` owns the
indexed sources in a `FixedIons` record. Result construction checks the
source indices and charge values, including when a failed or cancelled result carries effective
provenance.
Failed and cancelled results must not contain charge assignments.

The status is one of `success`, `invalid_input_or_request`, `no_executable_plan`, `numerical_failure`, or
`cancelled`. Invalid native API inputs generally throw `std::invalid_argument`; the owned facade converts
calculation failures and cooperative cancellation into result statuses. The lower-level
`CalculationRequest` overload propagates calculation and cancellation exceptions.

A `ChargeSet` owns the method ID, optional parameter-set ID, and source-indexed assignments.
Geometry-dependent methods return one assignment per molecule conformer. Geometry-independent methods
return one assignment per molecule. Every `ChargeAssignment` carries a molecule index and an optional
conformer index; its charge vector follows source atom order.

## Reuse, progress, and concurrency

`calculate(assessment, plan)` executes a concrete plan without repeating preparation, classification, or
applicability checks. Assessments and plans can be reused, and independent calculations can run
concurrently.

Derive from `calculation::CalculationObserver` to receive computation, target, and reduced-fragment
progress or to request cooperative cancellation. The observer is non-owning and callbacks may run on
oneTBB workers, so implementations must be thread-safe and must not throw. No partial charge set is
returned on cancellation.

## Molecular adapters

Public headers are provided under `chargefw/adapters`, with format-specific APIs in its `native` and
`gemmi` namespaces. The `all.h` headers are convenience includes for applications.

- `MolReader`, `SdfReader`, and `Mol2Reader` import native molecular formats.
- `JsonReader` imports ChargeFW molecule JSON 1.0.
- `PdbReader` and `MmcifReader` import Gemmi structures with explicit record, bond, and conformer policy.
- `native::mol2_output::Mol2Writer` emits generated MOL2 from an owned calculation result.
- `gemmi::mmcif_output::MmcifWriter` emits a fresh mmCIF document from an owned calculation result.
- `JsonWriter` serializes `ChargeCalculationResult` directly with the fixed `ChargeFW` generator identity;
  application-specific execution metrics are optional. `PortableId` retains absent, string, or signed
  64-bit integer record IDs without implicit stringification.

`RequestedCalculationProvenance` defaults to `execution_kind = "auto"`, matching automatic execution
selection. Set it explicitly when assembling a result for a different requested execution mode.

Readers return `ImportedMoleculeRecord`, which keeps the molecule together with source identity, import
diagnostics, and optional `MoleculeImportMetadata`. Its `SourceAtomReference`,
`SourceConformerReference`, `SourceStructuralLabels`, and `SourceConnectivity` values own their source
tokens and remain valid after reader destruction. The [molecular format reference](../reference/formats.md) documents
the populated mappings, supported subsets, reader policy, and generated output independently of the C++
types. `MoleculeImportMetadata::components` contains `SourceComponentInstance` values with a canonical
component ID and calculation-atom membership. Reader cursors are movable but not copyable; a moved reader
still requires its input stream to remain alive until reading is complete.
