# Calculate dipole moments

This script calculates charges for every conformer of a molecule with coordinates and prints each
conformer's net charge and dipole moment of the point charges, in debye:

```console
$ python calculate_dipole_moments.py paracetamol.sdf --format sdf --method eem
Model: eem/EEM_Ionescu2013_mpa_gas
molecule  conformer  net charge   |mu| (D)     mu_x      mu_y      mu_z
       0          0      +0.000      5.796    -4.969    +2.548    +1.551
```

With a method that does not use geometry, every conformer has the same charges, and only the geometry
changes the dipole. The script's `dipole_moment()` function takes charges, coordinates, and atomic
numbers, and can be used on its own.

Positions are measured from the center of nuclear charge. The dipole of a neutral molecule does not
depend on that origin, but the dipole of a charged molecule does, so compare charged molecules only with
the same origin. The dipole reflects the charge model: models fitted to population-analysis charges can
differ substantially from experimental dipoles or from charges fitted to the electrostatic potential.
Compare models and conformers with each other rather than treating one value as a measurement.

Download {download}`calculate_dipole_moments.py`.

```{literalinclude} calculate_dipole_moments.py
:language: python
```
