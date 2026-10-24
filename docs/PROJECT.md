# Project design

## Purpose

ChargeFW calculates empirical partial atomic charges from an explicit molecular graph, formal charges,
and, when required by a method, Cartesian coordinates. The same native engine is used by the C++ API,
the CLI, and the Python package.

ChargeFW is library-first. File readers and toolkit adapters translate external representations into an
owned, toolkit-neutral molecule model; method applicability, parameter matching, execution selection,
and calculation remain in the native library.

ChargeFW implements finite, non-periodic molecular calculations. Periodic and Ewald variants are not
implemented.

## Design principles

- **Explicit molecular input.** ChargeFW does not silently change topology, protonation, formal charge,
  hydrogen count, coordinates, method options, parameter sets, or execution mode.
- **Scientific applicability is separate from execution policy.** Missing coordinates, topology,
  element properties, or parameter coverage are hard applicability failures. Resource thresholds only
  guide automatic execution.
- **Selection is deterministic.** Unless constrained by the caller, methods and parameter sets are
  ordered by priority and then by stable ID.
- **Full and reduced execution are distinct.** Full execution follows the implemented finite-molecule
  method. Cutoff and cover are explicit approximations and are recorded in result provenance.
- **Mappings survive every boundary.** Molecule, atom, and conformer order are preserved in calculation
  results and adapters retain available source identities.
- **Algorithms are reusable and concurrent.** Prepared assessments and concrete plans own the state they
  need and can be executed repeatedly without global mutable calculation state.

## Architecture

```text
files, C++ values, or Python arrays
                 |
                 v
       toolkit-neutral molecules
       atoms + bonds + conformers
                 |
                 v
       prepared topology/geometry
                 |
        +--------+--------+
        |                 |
        v                 v
 bundled parameters   method registry
        |                 |
        +--------+--------+
                 v
       applicability assessment
                 |
                 v
       concrete execution plans
                 |
                 v
          full/cutoff/cover
                 |
                 v
      source-ordered charge sets
```

The public native namespaces mirror these responsibilities:

| Namespace | Responsibility |
| --- | --- |
| `chargefw::core` | Owned atoms, bonds, molecules, conformers, and periodic-table data |
| `chargefw::features` | Prepared topology, geometry, and spatial fragments |
| `chargefw::parameters` | Parameter data, loading, and atom/bond classification |
| `chargefw::methods` | Method metadata, options, requirements, applicability, and algorithms |
| `chargefw::calculation` | Assessment, deterministic planning, execution policy, progress, and calculation |
| `chargefw::charges` | Source-indexed charge assignments and charge sets |
| `chargefw::adapters` | Native and Gemmi-backed molecular input and output |

Application-facing calculation functions compose the owned assessment facade rather than exposing the
lower-level prepared-feature and classification interfaces. Callers need an explicit assessment only when
they want to inspect or reuse its concrete plans.

## Methods and parameters

The built-in registry contains 22 method IDs:

```text
abeem, charge2, delre, denr, dummy, eem, eqeq, eqeqc, formal, gdac,
kcm, mgc, mpeoe, peoe, qeq, sfkeem, smpqeq, sqe, sqeq0, sqeqp, tsef, veem
```

Bundled JSON parameter sets are installed with the library. Parameter-dependent methods are assessed
against this catalog. The [parameter-set JSON reference](PARAMETERS.md) documents its data format and
classifier semantics.

Methods declare their coordinate, topology, element-property, formal-charge, and parameter
requirements. Assessment reports both runnable plans and structured reasons for rejected candidates.
Method options are validated against each method's schema before execution.
Method metadata may also provide a concise human-readable note for important implementation scope or
input behavior. Notes explain structured requirements and deliberate variants; applicability checks remain
authoritative for a specific molecule and parameter set.

Method resource requirements describe asymptotic time and memory scaling in Big-O notation, where `n`
is the atom count and `m` is the bond count. Each term keeps the leading-order growth of one target's
full calculation; memory terms count the method's working storage beyond its inputs and output charges.
These terms guide execution planning; they do not predict a specific runtime or memory allocation.

Each method also estimates the cost of one target of a given atom and bond count: its working memory in
bytes, counted as for the memory term, and its operation count. Dense linear algebra is counted in
floating-point operations; other loops are counted per atom, bond, or atom-pair visit. The estimates
follow the dense allocations and loops of the implementation, so their leading terms match the declared
complexity. They are approximate: iterative methods are costed at their configured or maximum iteration
count, and operation counts compare work rather than predict seconds. Plans combine these estimates for
a whole calculation, as described under [assessment and execution](#assessment-and-execution).

PEOE and MPEOE initialize charge transfer from zero by default, so in that mode every connected
component must be formally neutral. Their `initial_charges=formal` option instead starts each atom at
its supplied formal charge and supports charged molecules while conserving the molecular formal-charge
total. Every connected component must be formally neutral for DelRe and ordinary SQE because their
implemented bond-transfer formulations start from zero and conserve zero charge within each component.
DENR, GDAC, KCM, and MGC have the same per-component requirement. Connected, net-neutral zwitterions
remain supported. SQE+q0 and SQE+qp include initial charges and support charged molecules when their
other requirements are met. SQE+q0 preserves the formal-charge seed total of each active connected
component; SQE+qp normalizes its fitted reference charges globally to the active target.

Charge2 implements the neutral-molecule inductive scheme from Abraham, Griffiths, and Loftus (1982).
Its one-, two-, and three-bond effects are applied as conservative bond transfers, and its iterative
polarizability feedback stops at the published charge tolerance or the configured maximum iteration count.
Every connected component must be formally neutral. The resonance and hybridization extensions from the
later Abraham and Hudson method are not implemented.

ABEEM places each bond charge at the point whose distances from the two nuclei are proportional to the
corresponding covalent radii, following Yang and Wang (1997). These centers determine nonincident
atom-bond and bond-bond interactions; the final atomic charges receive half of each incident bond charge.

QEq uses fixed, eV-valued atomic parameters in a single constrained solve with empirical two-center
Coulomb terms evaluated from angstrom coordinates. It does not iterate charge-dependent hydrogen. The
default DasGupta-Huzinaga term follows Oda and Hirono (2003); the optional Louwen-Vogt term uses the
geometric-mean approximation of Njo, Fan, and van de Graaf (1998).

Direct EEM calculations may include fixed Cartesian point sources. Their potential at active atom `i`
is `phi_i = sum_F(kappa * q_F / r_iF)`, using the selected parameter set's common `kappa`, so the
constrained-system right-hand side is `-A_i - phi_i`. The target charge remains the caller-supplied
active-atom total. Sources remain fixed while active charges respond; source polarization and charge
transfer between the active and fixed regions are excluded. Source/active distances must be valid and
nonzero. This implemented coupling does not establish quantitative accuracy for coordinated ions.

Full-mode SQE, SQE+q0, and SQE+qp accept fixed point sources. At active atom `i`, each source produces
the potential `phi_iF = q_F * erf(r_iF / (sqrt(2) * abs(width_i))) / r_iF`, with `q_F / r_iF` at zero
width. The potential is subtracted from the atom right-hand side before bond-transfer projection, with
no extra Coulomb factor. For SQE+qp, the parameterized `q0` reference charges are normalized once
globally to the active target charge; the source field changes active transfers without changing that
budget.

## Assessment and execution

Assessment prepares a molecule collection once, performs parameter classification, evaluates scientific
requirements, and expands applicable candidates into concrete plans. A plan contains the selected
method, parameter set, validated options, execution policy, and any policy warnings.

The CLI and Python facade opt into fixed ions with named imported components. The bundled catalog
contains 80 released monatomic CCD ions, with template-authoritative formal charges:

| Charge | Component IDs |
| --- | --- |
| -1 | `BR`, `CL`, `F`, `IOD` |
| +1 | `AG`, `AU`, `CS`, `CU1`, `D8U`, `K`, `LI`, `NA`, `RB`, `RH`, `TL` |
| +2 | `0BE`, `BA`, `CA`, `CD`, `CO`, `CU`, `EU`, `FE2`, `HG`, `MG`, `MN`, `NI`, `PB`, `PD`, `PT`, `RHF`, `SR`, `Y1`, `YB2`, `ZN` |
| +3 | `3CO`, `3NI`, `AL`, `AM`, `AU3`, `BS3`, `CE`, `CF`, `CR`, `CU3`, `DY`, `ER3`, `EU3`, `FE`, `GA`, `GD3`, `HO3`, `IN`, `IR3`, `LA`, `LU`, `MN3`, `ND`, `OS`, `PR`, `RH3`, `RU`, `SB`, `SM`, `TB`, `V`, `YB`, `YT3`, `ZCM`, `ZTM` |
| +4 | `4MO`, `4PU`, `4TI`, `IR`, `OS4`, `PT4`, `TH`, `ZR` |
| +6 | `6MO`, `W` |

The common-ion subset is `NA`, `K`, `MG`, `CA`, `CL`, `ZN`, `FE`, and `FE2`. Python exposes that subset
and the complete catalog as [immutable name selections](PYTHON.md#calculation-policy).
Each selected component must contain exactly one atom of the catalog element and have no incident
graph bond. Atom names do not affect recognition: component `CA` with one calcium atom supplies +2.
A selected template supplies its charge even when the imported atom has
a different formal charge. Input atoms remain unchanged. These selections are explicit caller intent,
independent of candidate parameter coverage; molecular ions such as sulfate remain part of the active graph.

Names are exact; repeated names are idempotent, and a known preset absent from the collection selects no
source. Unknown names and matching components that fail identity or graph-isolation validation reject the
request. Structural component identity
and canonical name selection follow [the format contract](FORMATS.md#pdb-and-mmcif-input). Each interface
resolves its named request to the indexed sources accepted by native `AssessmentRequest`; methods receive the
corresponding source positions and charges. See the [native API](NATIVE.md#assessment).

Assessment owns both the original collection and an active-molecule partition, then classifies capable methods
against that active collection. Rejection indices map back to the original molecules. EEM, SQE, SQE+q0, and
SQE+qp can be planned and executed in full, cutoff, or cover mode when their scientific prerequisites and
active charge budgets are satisfied. Reduced fragments are built only from active atoms; every fragment
receives the complete prescribed source field, including
sources outside its radius. Fixed atoms never become fragment centers, pivots, or charge variables.
Conservation corrections apply to active results, then the facade restores source charges exactly once in
original atom order. Effective provenance retains the resolved source indices and values;
the effective method, parameter set, and options define
the coupling.
This opt-in point-source approximation does not validate the response of coordinated metals.
Automatic selection considers warnings for each supported mode and may choose a reduced plan when full
execution carries a resource warning. Explicit selection can execute with its warning. Existing resource
thresholds continue to use active molecule size per target.

The active total is the formal-charge sum of unselected atoms in the prepared active molecule. Imported
formal charges on selected source atoms do not enter it; fractional source values contribute only at their
fixed sites. The modeled molecular total is `Q_model = Q_active + sum(fixed source charges)` and may
differ from the input formal-charge sum.
Assessment and execution do not change input atoms or infer chemical preparation. EEM constrains its solve
to the active total, and SQE+qp globally normalizes its fitted seeds to that target; SQE requires its
existing neutral active connected components, and SQE+q0 retains the active formal-charge seeds and
component totals.

The effective execution mode is always one of:

- `full`: calculate the complete molecular target;
- `cutoff`: solve one radius fragment per active atom and retain its mapped center charge;
- `cover`: solve radius fragments around active-atom-order pivots and retain charges in covered interiors.

Automatic planning prefers full execution. For methods classified as expensive, the default resource
policy moves to cutoff above 20,000 atoms and to cover above 80,000 atoms when those modes are supported.
Automatic reduced execution uses a 12 Å radius. Explicit reduced execution requires a radius of at least
8 Å. This lower bound was selected empirically because smaller fragments showed large RMS deviations from
full calculations. Larger radii generally improve agreement with full execution by including more of the
molecular environment, but increase calculation time. Explicit full execution may override a resource
threshold and returns a warning rather than changing the requested mode.

A thread limit caps the parallelism oneTBB allows the process, which is the hardware concurrency unless
the embedding application restricts it. Zero, or any larger limit, uses that full parallelism.

A plan can estimate its peak working memory and total operation count before execution, for a thread
limit. Full execution runs whole targets in parallel, so the memory estimate adds the workspaces of the
largest targets that can run together. Reduced execution runs targets one at a time and their fragments
in parallel: cutoff runs one fragment per active atom and cover one per pivot. Fragment sizes follow an
atom density of 0.1 atoms/Å³, typical of packed biomolecular structures with hydrogens, and never exceed
their target. Both modes add the retained output charges. The estimate describes calculation working
memory on top of the prepared assessment, not total process memory, and is most representative for
compact structures.

For ordinary calculations, cutoff and cover are implemented for:

```text
abeem, eem, eqeq, eqeqc, qeq, sfkeem, sqe, sqeq0, sqeqp
```

Before spatial cutting, reduced execution assigns every active atom a method-specific reference charge.
Each fragment inherits the sum of its active atoms' references as its target charge. After fragment
results have been assembled, a uniform correction restores the reference total separately in each group
whose charge the full method preserves:

| Methods | Source atom reference | Conserved group |
| --- | --- | --- |
| ABEEM, EEM, EQeq, EQeq+C, QEq, SFKEEM | Prepared active formal-charge total divided by active atom count | Entire calculation target |
| SQE | Zero | Each connected component in the pre-fragment active bond graph |
| SQE+q0 | Atomic formal charge | Each connected component in the pre-fragment active bond graph |
| SQE+qp | Parameterized charge, normalized once to the prepared active target's formal-charge total | Each connected component in the pre-fragment active bond graph |

The component labels always come from the pre-fragment active graph, not from pieces disconnected by an
individual spatial cut. Components remain together in a fragment when they are spatially nearby, allowing their
electrostatic interactions to alter atomic charges without changing the component reference totals.
SQE+qp retains its global normalization convention, so a component's preserved reference total need not
equal that component's formal charge.

Reduced calculations remain approximations. Cutting a bond removes transfers through the omitted region,
and omitted atoms no longer contribute their interactions to a fragment. Restoring conserved totals does
not reconstruct either effect, and the empirical radius guidance is not a general accuracy guarantee
across methods and molecular systems.

## Results and provenance

Geometry-dependent methods produce one assignment for each molecule conformer. Geometry-independent
methods produce one assignment per molecule. Every assignment identifies its source molecule and, when
applicable, conformer; atom values remain in source order.

Application-facing results distinguish success, invalid input or request, no executable plan, numerical
failure, and cancellation. Successful results retain the effective method, parameter set, validated
options, execution mode, radius, and execution warnings. The CLI additionally records requested policy,
diagnostics, timings, timestamps, and peak resident memory in its JSON result.

Cancellation is cooperative. A request is checked before each molecule or conformer target and before
each cutoff or cover fragment, so reduced calculations stop promptly, while a full calculation of one
target completes its solve before cancellation takes effect. A cancelled result contains no charges.

## Molecular data scope

The language-independent [molecular format reference](FORMATS.md) defines supported inputs and outputs,
exactly what each reader imports, and how each writer preserves or generates molecular data.

Imported records keep adapter metadata beside the normalized molecule rather than inside
`core::Molecule`. The metadata owns only source correspondence, record-local import policy, diagnostics,
and a small connectivity summary; it does not retain source documents, parser objects, normalized-value
snapshots, or per-bond provenance.

The application-facing charge result owns those ordered imported records together with the toolkit-neutral
execution result. Output adapters consume that boundary directly, keeping source correspondence out of the
scientific calculation engine without requiring parallel caller-supplied structures or contexts.

ChargeFW is not a general chemistry-preparation toolkit. It does not provide SMILES parsing, arbitrary
bond perception, sanitization, protonation, hydrogen addition or removal, conformer generation, or
geometry optimization. Structural template bonding covers common amino acids, nucleotides, water, and
basic polymer links; it is not a complete Chemical Component Dictionary provider.

## Versions

The native library, the command-line application, and the Python package share one
[PEP 440](https://peps.python.org/pep-0440/) version, which continues ChargeFW2's lineage. Pre-releases
such as `3.0.0a1` may change any interface. From 3.0.0, incompatible changes to documented interfaces
wait for a new major version, and every change that affects calculated charges is listed in the
[changelog](https://github.com/sb-ncbr/chargefw/blob/main/CHANGELOG.md). New parameter sets can change
automatic selection, so name the method and parameter set explicitly when results must stay the same
across releases.
