# Large structures

Most geometry-dependent methods solve a system of linear equations for the whole molecule. Their time
grows with the cube of the number of atoms and their memory with its square, so a full calculation of a
large protein complex can take hours or not fit into memory. For such structures, ChargeFW can calculate
charges from overlapping fragments instead, with two approximations called cutoff and cover.

## What ChargeFW does by default

ChargeFW calculates the whole molecule, called full execution, unless the molecule is large and the
method expensive. A method is expensive when its time grows with the cube of the atom or bond count, or
its memory with the square; the **Time** column of the [method table](methods.md#available-methods)
shows the time. For expensive methods, ChargeFW chooses:

| Atoms in the molecule | Execution |
| --- | --- |
| up to 20,000 | full |
| 20,001 to 80,000 | cutoff, with a 12 Å radius |
| more than 80,000 | cover, with a 12 Å radius |

Methods that are not expensive, such as PEOE, always run in full. Expensive methods without cutoff and
cover, such as DelRe or SMP/QEq, are not chosen automatically for molecules above 20,000 atoms; request
full execution to use them anyway. The result records the execution mode and radius that were used.

## Choosing the execution yourself

Request a mode and, for cutoff and cover, a radius of at least 8 Å:

```bash
chargefw calculate --execution cover --radius 14 structure.cif output
```

```python
result = chargefw.calculate(molecules, method="eem", execution="cover", radius=14.0)
```

With automatic execution, a radius sets the radius ChargeFW uses if it switches to cutoff or cover.

Explicit full execution also runs above the size limits; the result then carries a warning. The limits
themselves can be changed with `--cutoff-atom-threshold` and `--cover-atom-threshold` on the command line,
or `cutoff_threshold` and `cover_threshold` in Python. `unlimited` on the command line, or `None` in
Python, removes a limit. The cover limit must not be lower than the cutoff limit, and removing the
cutoff limit requires removing the cover limit too.

Before a long calculation, `chargefw applicability` estimates the working memory and the number of
operations of every possible plan. For crambin with hydrogens (642 atoms) and SQE+qp:

```text
plan method=sqeqp parameter_set=SQEqp_Schindler2021_CCD_gen execution=full memory=6.4MiB operations=1.9e+08
plan method=sqeqp parameter_set=SQEqp_Schindler2021_CCD_gen execution=cutoff memory=102.2MiB operations=1.2e+11
plan method=sqeqp parameter_set=SQEqp_Schindler2021_CCD_gen execution=cover memory=102.2MiB operations=1.1e+10
```

For a molecule this small, full execution is by far the cheapest. Cover needs about ten times fewer
operations than cutoff, because it calculates far fewer fragments.

## How cutoff and cover work

Both approximations were introduced as EEM Cutoff and EEM Cover in AtomicChargeCalculator
([Ionescu 2015](https://doi.org/10.1186/s13321-015-0099-x)). ChargeFW uses them for more methods, and
its cover selects fragments and assigns charges differently, as described below.

**Cutoff** builds one fragment for every atom: the atom and all atoms within the radius around it. It
calculates the fragment and keeps the charge of the central atom only.

**Cover** goes through the atoms in input order. Each atom that no fragment covers yet becomes the
center of a new fragment and covers all atoms within 3 Å of it that are not covered yet. The fragment
contains all atoms within the radius around the center, and the covered atoms keep the charges
calculated for them in this fragment. In AtomicChargeCalculator, the fragment centers were instead
chosen so that no two centers are bonded and every atom has a center within two bonds, and the charge of
an atom was summed from every fragment that contained it.

Before cutting the molecule, ChargeFW gives every atom a reference charge: its share of the molecule's
charge, matching what the method keeps constant in a full calculation. Each fragment receives the sum of
the reference charges of its atoms as its total charge. After all fragments are calculated, ChargeFW adds
the same correction to every atom of a group, so that the group's total equals the sum of its reference
charges:

| Methods | Reference charge of an atom | Group whose total is kept |
| --- | --- | --- |
| EEM, QEq, ABEEM, SFKEEM, EQeq, EQeq+C | Total formal charge divided by the number of atoms | The whole molecule |
| SQE | Zero | Each bonded part of the molecule |
| SQE+q0 | Its formal charge | Each bonded part of the molecule |
| SQE+qp | Its parameterized charge, shifted so that all of them add up to the total formal charge | Each bonded part of the molecule |

For SQE+qp, the total kept for a bonded part can therefore differ from the part's formal charge, as in a
full SQE+qp calculation. The size of the correction depends on how well the fragments represent the
molecule; it is not guaranteed to be small.

Both approximations leave out what lies beyond the radius: bonds cut at the edge of a fragment and atoms
farther away. [Fixed ions](fixed-ions.md) are the exception: they never become fragment centers, every
fragment includes all of them, also those beyond the radius, and their charges are not corrected. A
larger radius brings the charges closer to a full calculation and makes the calculation slower.
Below 8 Å, the charges differed too much from full calculations, which is why 8 Å is the smallest
radius allowed. To judge the error for your system, calculate a smaller but similar structure both
in full and with the approximation, and compare.

## Threads and estimates

`--threads` on the command line, or `threads` in Python, limits the number of threads. By default,
ChargeFW uses all available processor cores. Full execution calculates several molecules in parallel;
cutoff and cover calculate one molecule at a time and its fragments in parallel.

The estimates count the working memory of the calculation itself, not the memory of the loaded
structure, and the operations of the method; they compare plans and do not predict seconds. Fragment
sizes are estimated for the atom density of packed proteins with hydrogens, about 0.1 atoms per Å³, so
the estimates are most accurate for compact structures.
