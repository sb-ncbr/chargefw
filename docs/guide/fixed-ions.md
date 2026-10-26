# Fixed ions

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
and the complete catalog as [immutable name selections](../interfaces/python.md#calculation-policy).
Each selected component must contain exactly one atom of the catalog element and have no incident
graph bond. Atom names do not affect recognition: component `CA` with one calcium atom supplies +2.
A selected template supplies its charge even when the imported atom has
a different formal charge. Input atoms remain unchanged. These selections are explicit caller intent,
independent of candidate parameter coverage; molecular ions such as sulfate remain part of the active graph.

Names are exact; repeated names are idempotent, and a known preset absent from the collection selects no
source. Unknown names and matching components that fail identity or graph-isolation validation reject the
request. Structural component identity
and canonical name selection follow [the format contract](../reference/formats.md#pdb-and-mmcif-input). Each interface
resolves its named request to the indexed sources accepted by native `AssessmentRequest`; methods receive the
corresponding source positions and charges. See the [native API](../interfaces/cpp.md#assessment).

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

## Point-source terms

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
