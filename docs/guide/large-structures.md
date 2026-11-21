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

| Atoms in the largest molecule | Execution |
| --- | --- |
| up to 20,000 | full |
| 20,001 to 80,000 | cutoff, with a 12 Å radius |
| more than 80,000 | cover, with a 12 Å radius |

Methods that are not expensive, such as PEOE, always run in full. Expensive methods without cutoff and
cover, such as DelRe or SMP/QEq, are not chosen automatically for molecules above 20,000 atoms; request
full execution to use them anyway. The result records the execution mode and radius that were used.

ChargeFW chooses one execution for all molecules of a calculation, so the largest molecule decides: a
50-atom molecule calculated together with a 30,000-atom protein gets cutoff too. To keep small molecules
in full execution, calculate them separately from large ones.

## How long a calculation takes

These times were measured with the command line on an 8-core desktop processor (AMD Ryzen 7 9700X,
8 threads, the default), including reading the input and writing the result:

| Structure | Atoms | Method | Execution | Time | Peak memory |
| --- | --- | --- | --- | --- | --- |
| [1L2Y](input.md#check-the-molecule-before-calculating) | 304 | EEM or SQE+qp | full | 0.01 s | 10 MB |
| 10AW, without hydrogens | 10,479 | EEM | full | 1.5 s | 0.9 GB |
| 10AW, without hydrogens | 10,479 | SQE+qp | full | 1.0 s | 1.8 GB |
| 8YAX | 40,737 | EEM | full | 69 s | 13 GB |
| 8YAX | 40,737 | SQE+qp | full | 40 s | 26 GB |
| 8YAX | 40,737 | EEM | cutoff, 12 Å (default) | 6.8 s | 0.10 GB |
| 8YAX | 40,737 | EEM | cover, 12 Å | 1.2 s | 0.10 GB |
| 8YAX | 40,737 | SQE+qp | cutoff, 12 Å (default) | 10.6 s | 0.14 GB |
| 8YAX | 40,737 | SQE+qp | cover, 12 Å | 1.7 s | 0.14 GB |

Four times as many atoms made the full calculation about 40 times slower and its memory 14 times
larger. Cutoff and cover grow about in proportion to the number of atoms. On this structure they were
4 to 60 times faster than the full calculation and needed less than a hundredth of its memory.

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
themselves can be changed with `cutoff_threshold` and `cover_threshold` in Python, or
`--cutoff-threshold` and `--cover-threshold` on the command line. `unlimited` removes a limit. The
cover limit must not be lower than the cutoff limit, and removing the cutoff limit requires removing
the cover limit too.

Before a long calculation, `chargefw applicability` estimates the memory and the number of
operations of every plan: each combination of method, parameter set, and execution mode that
applies. For the Trp-cage miniprotein {download}`1l2y-h.cif` (304 atoms) and SQE+qp:

```console
$ chargefw applicability --method sqeqp --parameter-set SQEqp_Schindler2021_CCD_gen 1l2y-h.cif
plan method=sqeqp parameter_set=SQEqp_Schindler2021_CCD_gen execution=full memory=1.4MiB operations=2.0e+07
plan method=sqeqp parameter_set=SQEqp_Schindler2021_CCD_gen execution=cutoff memory=23.0MiB operations=6.1e+09
plan method=sqeqp parameter_set=SQEqp_Schindler2021_CCD_gen execution=cover memory=23.0MiB operations=5.7e+08
...
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
| SQE+qp | Its parameterized charge, shifted so that each bonded part adds up to its formal charge | Each bonded part of the molecule |

Each fragment of an SQE+qp calculation starts from these reference charges, as the full calculation
does. The size of the correction depends on how well the fragments represent the molecule; it is not
guaranteed to be small.

Both approximations leave out what lies beyond the radius: bonds cut at the edge of a fragment and atoms
farther away. [Fixed ions](fixed-ions.md) are the exception: they never become fragment centers, every
fragment includes all of them, also those beyond the radius, and their charges are not corrected. A
larger radius brings the charges closer to a full calculation and makes the calculation slower.
Below 8 Å, the charges differed too much from full calculations, which is why 8 Å is the smallest
radius allowed. To judge the error for your system, calculate a smaller but similar structure both
in full and with the approximation, and compare.

## Threads and estimates

`--threads` on the command line, or `threads` in Python, limits the number of threads. By default,
ChargeFW uses one thread per processor core that it is allowed to use. A core that runs two threads at
once (hyper-threading) counts once: a second thread per core made cutoff and cover 10 to 20% faster in
our measurements, and you can ask for it with a higher number.

Full execution calculates several molecules in parallel. An expensive method also uses the threads
within one large molecule: when its system of equations has 1,500 or more unknowns, which for most
methods means about 1,500 atoms. Cutoff and cover calculate one molecule at a time and its
fragments in parallel. The charges do not depend on the number of threads.

More threads shorten the full calculation of one molecule, but each added thread helps less. For 10AW
(10,479 atoms) on the processor above:

| Threads | EEM | SQE+qp |
| --- | --- | --- |
| 1 | 6.0 s | 3.8 s |
| 2 | 3.2 s | 2.1 s |
| 4 | 2.0 s | 1.4 s |
| 8 | 1.5 s | 1.0 s |

### On a shared server

When several calculations run on one machine, set `threads` for each of them instead of relying on the
default, which would give every calculation all cores:

- Choose the thread count so that the number of calculations running at once, times their threads,
  is at most the number of processor cores.
- Prefer more calculations with 2 to 4 threads each over fewer with many. In the table above, 2 threads
  are almost twice as fast as one and 4 threads about three times, while 8 threads are four times
  faster.
- Check the memory before you start a calculation. The [estimate](#estimates) of a plan gives it
  without calculating. Ask for it with the thread count you will calculate with, such as
  `plan.estimate(threads=4)` in Python or `chargefw applicability --threads 4` on the command line.
  Memory depends on the method and grows with the square of the atom count. As examples, the full
  calculation of a 20,000-atom molecule needs about 3 GB with EEM and 7 GB with SQE+qp; other methods can
  need twice that.
- Leave the execution automatic for requests you do not control. A request for full execution of a
  40,000-atom structure took 13 GB with EEM and 26 GB with SQE+qp in the
  [table above](#how-long-a-calculation-takes).
- To limit a container, give it a processor set, such as `--cpuset-cpus` in Docker. ChargeFW then
  counts only the cores of that set.

### Estimates

The estimates count the memory of the calculation itself, not that of the loaded structure, and the
operations of the method; they compare plans and do not predict seconds. Fragment sizes are
estimated for the atom density of packed proteins with hydrogens, about 0.1 atoms per Å³, so the
estimates are most accurate for compact structures.