# Calculate with a fixed ion

This script calculates a PDB or mmCIF structure with selected ions kept as
[fixed charges](../guide/fixed-ions.md), and reports the total charge and the charges around each ion.
For [water-ions.cif](../guide/water-ions.cif):

```console
$ python calculate_with_fixed_ions.py water-ions.cif --format mmcif --ion MG --ion CA
molecule=0 conformer=0: total=+4.000 e (active atoms +0.000, fixed ions +4.000)
  fixed I/MG 9/MG charge=+2.0 e
    W/HOH 5/O            O   distance=3.00 A charge=-0.877 e
    I/CA 10/CA           Ca  distance=3.00 A charge=+2.000 e
  fixed I/CA 10/CA charge=+2.0 e
    I/MG 9/MG            Mg  distance=3.00 A charge=+2.000 e
```

Each fixed ion is listed by chain, residue, and atom name, followed by the atoms within `--radius`
(default 3.0 Å), such as the atoms coordinating it. Select ions with repeated `--ion`, or with
`--common-ions` for the common ions; ions that are not in the structure are ignored. `--result-json`
also writes the [result JSON](../guide/results.md).

The script uses SQE+qp with `SQEqp_Schindler2021_CCD_gen` unless you choose another method or
parameter set; only methods that support fixed ions can be used.

Download {download}`calculate_with_fixed_ions.py`.

```{literalinclude} calculate_with_fixed_ions.py
:language: python
```
