# C++ library

ChargeFW is a C++23 shared library for Linux, used through the `chargefw::core` CMake target. An
assessment keeps its own copies of the molecules and parameter sets, and results hold their own charges.
The one lifetime to manage is that of an assessment: its plans refer to it, as described in
[Assessment and calculation](#assessment-and-calculation).

## Build and link

ChargeFW requires CMake 3.28 or newer, Ninja, a GCC or Clang toolchain with C++23 support, and installed
development packages for nlohmann/json 3.11.3 or newer, Eigen 3.4 or 5, nanoflann 1.7 or newer,
oneTBB 2022.1 or newer, Gemmi 0.7.4 or newer, and, for the command line, CLI11 2.4 or newer. To build
without them, set `CHARGEFW_FETCH_DEPENDENCIES=ON`; ChargeFW then downloads pinned, checksum-verified
sources of all of them.

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

Set `CHARGEFW_ENABLE_NATIVE_OPTIMIZATIONS=ON` only for an installation that stays on the build machine;
interprocedural optimization alone does not make the binaries CPU-specific.

Use the installed package from another CMake project:

```cmake
find_package(chargefw CONFIG REQUIRED)
target_link_libraries(my_application PRIVATE chargefw::core)
```

The installed package requires the exact nlohmann/json and Gemmi versions ChargeFW was built against,
because their types appear in the public headers. With installed dependency packages, your project needs
those packages discoverable and the oneTBB runtime library available. With
`CHARGEFW_FETCH_DEPENDENCIES=ON`, the installation also contains the nlohmann/json and Gemmi development
packages and the oneTBB runtime.

The installation can be moved as a whole: `parameters::load_default_parameter_sets()` finds the bundled
parameter files relative to the loaded library.

## Calculating charges

Build a molecule, assess it, and calculate:

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

`assess()` copies the molecules and parameter sets of an lvalue request and moves those of an rvalue
request.

## Reading files and writing results

The readers return `adapters::ImportedMoleculeRecord` values, which keep each molecule together with
where it came from. The writers take the records and the calculation together, as an
`adapters::ChargeCalculationResult`:

```cpp
#include <chargefw/adapters/charge_result.h>
#include <chargefw/adapters/native/json_output.h>
#include <chargefw/adapters/native/mol2_output.h>
#include <chargefw/adapters/native/sdf_input.h>
#include <chargefw/calculation/calculation.h>
#include <chargefw/config.h>
#include <chargefw/parameters/io/parameter_set_io.h>

#include <fstream>
#include <iostream>
#include <utility>
#include <vector>

int main() {
    using namespace chargefw;

    auto input = std::ifstream{"paracetamol.sdf"};
    auto reader = adapters::native::sdf_input::SdfReader{input, "paracetamol.sdf"};
    auto records = std::vector<adapters::ImportedMoleculeRecord>{};
    auto molecules = std::vector<core::Molecule>{};
    while (auto record = reader.next()) {
        molecules.push_back(record->molecule);
        records.push_back(std::move(*record));
    }

    auto request = calculation::AssessmentRequest{
        .molecules = core::MoleculeCollection{std::move(molecules)},
        .parameter_sets = parameters::load_default_parameter_sets(),
        .method_id = "eem",
        .parameter_set_id = "EEM_Cheminf_b3lyp_npa",
    };
    auto assessment = calculation::assess(std::move(request));
    auto execution = calculation::calculate(assessment);
    if (!execution.calculated()) {
        std::cerr << execution.failure_message.value_or("calculation failed") << '\n';
        return 1;
    }

    const auto requested = adapters::RequestedCalculationProvenance{
        .method_id = "eem", .parameter_set_id = "EEM_Cheminf_b3lyp_npa"};
    const auto result =
        adapters::make_charge_calculation_result(std::move(records), requested, std::move(execution));

    auto mol2 = std::ofstream{"paracetamol.mol2"};
    adapters::native::mol2_output::Mol2Writer{mol2}.write(result, CHARGEFW_VERSION_STRING);
    auto json = std::ofstream{"paracetamol.json"};
    adapters::native::json_output::JsonWriter{json}.write(result, CHARGEFW_VERSION_STRING);
}
```

| Header | Class |
| --- | --- |
| `chargefw/adapters/native/mol_input.h`, `sdf_input.h`, `mol2_input.h` | `MolReader`, `SdfReader`, `Mol2Reader` |
| `chargefw/adapters/native/json_input.h` | `JsonReader`, for molecule JSON 1.0 |
| `chargefw/adapters/gemmi/pdb_input.h`, `mmcif_input.h` | `PdbReader`, `MmcifReader`, with selection, bond, and conformer options |
| `chargefw/adapters/native/mol2_output.h` | `Mol2Writer` |
| `chargefw/adapters/native/pqr_output.h` | `PqrWriter`, for one structure read from PDB or mmCIF |
| `chargefw/adapters/gemmi/mmcif_output.h` | `MmcifWriter`, which writes a new mmCIF document with one or more charge sets |
| `chargefw/adapters/native/json_output.h` | `JsonWriter`, for result JSON; execution metrics are optional |

`MmcifWriter` also writes several charge sets for the same records into one file. Pass pointers to
results calculated from equal records, for example from copies of `records` above:

```cpp
auto cif = std::ofstream{"paracetamol.cif"};
adapters::gemmi::mmcif_output::MmcifWriter{cif}.write(std::array{&eem_result, &sqeqp_result},
                                                      CHARGEFW_VERSION_STRING);
```

The [mmCIF reference](../reference/formats.md#mmcif) describes when results are refused.

The writers record ChargeFW and the version you pass as the program that wrote the file; pass
`CHARGEFW_VERSION_STRING` from `<chargefw/config.h>`.

`chargefw/adapters/native/all.h` and `gemmi/all.h` include all headers of a namespace, for applications.
A reader can be moved but not copied, and its input stream must stay open until reading is complete. The
records stay valid after the reader is destroyed.

`ImportedMoleculeRecord` holds the molecule, its source identity, import diagnostics, and optional
`MoleculeImportMetadata`: the source atoms (`SourceAtomReference`), conformers
(`SourceConformerReference`), structural labels, connectivity, and residues
(`SourceComponentInstance`). The [format reference](../reference/formats.md) describes what each reader
fills in. `PortableId` holds a record or atom ID that is absent, a string, or a 64-bit integer.

`make_charge_calculation_result()` checks that the records fit the calculation: their number, their atom
and conformer counts, the assignments, and the method and parameter set. It cannot check that they are
the records you calculated, so pass the same records whose molecules you calculated.
`RequestedCalculationProvenance` records what was requested, for the result JSON; its `execution_kind`
defaults to `"auto"`, so set it when you request another execution mode.

## Namespaces

| Namespace | Contents |
| --- | --- |
| `chargefw::core` | Atoms, bonds, molecules, conformers, and the periodic table |
| `chargefw::features` | Prepared topology, geometry, and spatial fragments |
| `chargefw::parameters` | Parameter sets, loading, and atom and bond classification |
| `chargefw::methods` | Method metadata, options, requirements, applicability, and algorithms |
| `chargefw::calculation` | Assessment, plans, execution policy, progress, and calculation |
| `chargefw::charges` | Charge assignments and charge sets |
| `chargefw::adapters` | File readers and writers |

## Molecules

`core::Molecule` holds:

- a vector of `Atom` values in source order, each with atomic number, formal charge, and optional name;
- single, double, or triple `Bond` values between atom indices;
- zero or more `Conformer` values, each with one position per atom in ångströms; and
- an optional name.

`core::MoleculeCollection` holds molecules in order, with an optional name. Atomic numbers 1 to 100 are
accepted; methods and parameter sets may support fewer elements. The constructors check atomic numbers,
bond indices and duplicates, and the number of positions per conformer, and throw
`std::invalid_argument` on errors. Coordinates that are not finite are kept, but methods that need
geometry then reject the molecule.

## Assessment and calculation

`calculation::AssessmentRequest` contains:

- the `MoleculeCollection` and the parameter sets;
- an optional method ID and parameter-set ID;
- options per method, as `methods::MethodOptions`;
- strict or permissive parameter matching;
- an `ExecutionSelection`: automatic, full, cutoff, or cover, with an optional radius;
- a `ResourcePolicy` with the automatic size limits, where `std::nullopt` removes a limit; and
- optional `FixedIons`, described [below](#fixed-ions).

Parameter-set IDs in one request must be unique. `calculation::assess()` returns an `AssessmentResult`
with the `molecules()`, the possible `plans()` in selection order, the `rejections()` with reasons, the
`default_plan()`, and the time taken. A plan has its method, parameter set, all option values, an
`ExecutionPolicy` with the mode and radius, and warnings.
[Choosing a method](../guide/methods.md#choosing-a-method) and
[Large structures](../guide/large-structures.md) explain selection and execution.

`calculation::calculate(assessment)` runs the default plan, and `calculate(assessment, plan)` a chosen
one, without preparing or checking the molecules again. A plan, and every copy of it, refers to data
of the assessment that created it, so keep the assessment alive while you use its plans, and use them
only with that assessment. An assessment can be moved into a new object, but not assigned. Both take a thread limit, where `0` uses all available cores, and an optional observer, described
[below](#progress-and-cancellation). An assessment can be calculated repeatedly, also from several
threads at once:

```cpp
auto assessment = calculation::assess(std::move(request));
auto serial = calculation::calculate(assessment, 1);
auto parallel = calculation::calculate(assessment, 4);
```

`calculation::estimate(assessment, plan, max_threads)` returns a `methods::ResourceEstimate` with the
plan's approximate peak working memory and number of operations; [Threads and
estimates](../guide/large-structures.md#threads-and-estimates) explains them.

## Results and errors

`calculation::ExecutionResult` has an `ExecutionStatus`, an optional `charges::ChargeSet`, the
rejections, the model that was used in `effective`, an optional failure message, and the time spent on
assessment and calculation. The status is `success`, `invalid_input_or_request`, `no_executable_plan`,
`numerical_failure`, or `cancelled`; failed and cancelled results contain no charges.

A `ChargeSet` has the method ID, the optional parameter-set ID, and the assignments. Methods that use
geometry give one assignment per conformer, the others one per molecule. Each `ChargeAssignment` has a
molecule index, an optional conformer index, and the charges in atom order.
[Results](../guide/results.md) explains assignments and atom order.

Invalid arguments throw `std::invalid_argument`, also from `assess()`. A valid assessment reports
methods that do not apply as rejections, and `default_plan()` is null when none applies; it has no
status. The `calculate()` overloads that take an assessment report failures and cancellation through
the `ExecutionStatus`; the lower-level `CalculationRequest` overload throws them as exceptions
instead.

## Fixed ions

`AssessmentRequest::fixed_ions` keeps chosen atoms as fixed point charges, as described in
[Fixed ions](../guide/fixed-ions.md). Its `sources` pair a charge with a molecule index and an atom index
in the collection. Unlike the command line and Python, C++ takes the charges directly, not ion names from
the catalog. Leave it empty for ordinary calculations.

Each source must refer to a different atom, have a finite charge, and refer to an atom without bonds.
Every affected molecule needs at least one atom that is not a source, and at least one conformer, with
finite coordinates throughout. No atom that is not a source may sit exactly at a source's position. Invalid sources throw `std::invalid_argument`. With
sources, only methods that support fixed ions are assessed, and the formal charges of the source atoms do
not count toward the charge of the calculated atoms. A successful result contains the fixed charges at
their atoms, and `ExecutionResult::effective` records the sources.

Code that calls a method directly, through `Method::calculate()`, passes fixed charges as a
`std::span<const methods::FixedPointSource>` at the end of `methods::CalculationInput`; each source has a
position and a charge. The input refers to the sources without copying them, so the source elements
and their storage must stay valid for every call that uses the input; the span object itself need
not. Such
code must check `MethodRequirements::supports_fixed_point_sources` itself and pass an empty span to other
methods, set the target charge to the charge of the calculated atoms, supply finite charges and
coordinates with no source at an atom's position, and assemble the results; the method calculates only
the atoms it is given.

## Methods and parameter sets

`methods::method_registry()` lists the built-in methods by descending priority, then ID, with their
metadata, notes, requirements, and option schemas. `methods::complexity_notation()` formats the declared
time and memory complexity, and `Method::estimate()` returns the approximate memory and operations of
one molecule of a given `TargetSize`. A method supports cutoff and cover when
`ResourceRequirements::supports_cutoff` and `supports_cover` are set; its `reduced_charge_policy`
says how the fragments' charges are kept, as described in
[How cutoff and cover work](../guide/large-structures.md#how-cutoff-and-cover-work). A method that
declares common parameters must also declare atom or bond parameters; make standalone constants method
options instead, or assessment refuses the method.

`parameters::load_default_parameter_sets()` loads the bundled parameter sets, and
`parameters::load_parameter_sets_json_directory()` the sets in another directory. Both group the
sets by method ID and order each method's sets by descending priority, then ID.
`parameters::load_parameter_set_json_file()` loads one file, and `load_parameter_set_json()` one set
from a `std::istream`. These functions are declared in
`<chargefw/parameters/io/parameter_set_io.h>`. The [parameter-set
format](../reference/parameter-format.md) describes the files.

`<chargefw/parameters/io/parameter_data_paths.h>` declares
`parameters::bundled_parameter_directory()`, the directory with the bundled sets of the installed
library, and `parameters::default_parameter_directories()`, every directory
`load_default_parameter_sets()` searches.

## Progress and cancellation

Derive from `calculation::CalculationObserver` and pass it as the last argument of `calculate()` to
receive progress for the computation, each molecule or conformer, and the fragments of cutoff and cover,
or to cancel. The calculation does not take ownership of the observer. Callbacks may run on several
worker threads at once, so they must be thread-safe and must not throw. A cancelled calculation returns
no charges.
