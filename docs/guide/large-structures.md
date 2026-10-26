# Large structures

**Full and reduced execution are distinct.** Full execution follows the implemented finite-molecule
method. Cutoff and cover are explicit approximations and are recorded in result provenance.

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

## Threads and resource estimates

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
