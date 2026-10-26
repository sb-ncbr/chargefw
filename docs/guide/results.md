# Results

**Mappings survive every boundary.** Molecule, atom, and conformer order are preserved in calculation
results and adapters retain available source identities.

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
