# Analyze conformer-dependent charges

This script generates conformers of a molecule with RDKit, calculates EEM charges for each, and lists
the atoms whose charges vary most between the conformers:

```console
$ pip install --pre "chargefw[rdkit]"
$ python analyze_rdkit_conformer_charges.py "CCO" --conformers 20
Method: eem/EEM_NEEMP_ccd2016_npa
Conformers: 20 equally weighted, unoptimized ETKDGv3 structures
Atom          mean         std         min         max
H6         0.214496    0.015047    0.188914    0.239530
C2        -0.018663    0.011457   -0.037413    0.000001
...
```

The conformers are generated with a fixed random seed and are not optimized, so the result shows how
sensitive the charges are to geometry within one set of conformers. It is not a thermodynamic average,
and atoms with varying charges are not necessarily chemically important. The analysis function in the
script also accepts conformers you prepared yourself.

Download {download}`analyze_rdkit_conformer_charges.py`.

```{literalinclude} analyze_rdkit_conformer_charges.py
:language: python
```
