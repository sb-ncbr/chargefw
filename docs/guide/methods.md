# Methods and parameter sets

ChargeFW implements finite, non-periodic molecular calculations. Periodic and Ewald variants are not
implemented.

**Scientific applicability is separate from execution policy.** Missing coordinates, topology,
element properties, or parameter coverage are hard applicability failures. Resource thresholds only
guide automatic execution.

**Selection is deterministic.** Unless constrained by the caller, methods and parameter sets are
ordered by priority and then by stable ID.

## Methods and parameters

The built-in registry contains 22 method IDs:

```text
abeem, charge2, delre, denr, dummy, eem, eqeq, eqeqc, formal, gdac,
kcm, mgc, mpeoe, peoe, qeq, sfkeem, smpqeq, sqe, sqeq0, sqeqp, tsef, veem
```

Bundled JSON parameter sets are installed with the library. Parameter-dependent methods are assessed
against this catalog. The [parameter-set JSON reference](../reference/parameter-format.md) documents its data format and
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
a whole calculation, as described under [resource estimates](large-structures.md#threads-and-resource-estimates).

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

## Assessment and selection

Assessment prepares a molecule collection once, performs parameter classification, evaluates scientific
requirements, and expands applicable candidates into concrete plans. A plan contains the selected
method, parameter set, validated options, execution policy, and any policy warnings.
