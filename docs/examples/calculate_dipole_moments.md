# Calculate dipole moments

`calculate_dipole_moments.py` reads every conformer of any supported input
with coordinates, calculates charges with one model, and prints each conformer's net charge and
point-charge dipole moment in debye. Charges from a geometry-independent method apply to every
conformer, so only the geometry changes the dipole. The reusable `dipole_moment()` function takes source-
ordered charges, coordinates, and atomic numbers.

```bash
python calculate_dipole_moments.py ensemble.cif --format mmcif --method eem
```

Positions are measured from the center of nuclear charge. A neutral molecule's dipole does not depend on
that origin, but a charged molecule's dipole does, so compare charged species only with the same origin.
A point-charge dipole reflects the charge model: empirical models fitted to population-analysis charges
can differ substantially from experimental or electrostatic-potential-derived dipoles, so compare models
and conformers rather than treating one value as a measurement.

Download {download}`calculate_dipole_moments.py`.

```{literalinclude} calculate_dipole_moments.py
:language: python
```
