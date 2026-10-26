# Calculate with a fixed ion

`calculate_with_fixed_ions.py` reads a PDB or mmCIF structure, fixes
selected monatomic ion components as point charges, and calculates the remaining active atoms. For each
conformer it reports the charge balance between active atoms and fixed ions, then every fixed ion by its
author chain, residue, and atom name with the charges of atoms within `--radius` (default 3.0 Å), such as
the ion's coordinating atoms. `--result-json` also writes the complete result.

```bash
python calculate_with_fixed_ions.py structure.cif \
    --format mmcif --ion MG --ion CA --result-json charges.json
python calculate_with_fixed_ions.py structure.cif --format mmcif --common-ions
```

Select ions with repeated `--ion` or with `--common-ions` for the common ion catalog; only matching
components present in the input become fixed sources. The example uses `sqeqp` with
`SQEqp_Schindler2021_CCD_gen` by default; override the method or parameter set when the active structure
requires another model. Only methods that support fixed point sources can be used.

Download {download}`calculate_with_fixed_ions.py`.

```{literalinclude} calculate_with_fixed_ions.py
:language: python
```
