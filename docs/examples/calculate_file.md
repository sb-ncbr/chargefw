# Calculate a molecular file

This script reads a molecular file in any supported format, calculates charges with the method you
choose, and prints one row per atom. For [paracetamol.sdf](../getting-started/paracetamol.sdf):

```console
$ python calculate_file.py paracetamol.sdf --format sdf --method eem
molecule  conformer  atom_index  atom_id  atom_name  element      charge
       0          0           0        0  -          C         -0.592069
       0          0           1        1  -          C          0.662721
       0          0           2        2  -          O         -0.719530
...
```

A `-` in the conformer column marks charges from a method that does not use geometry. `--csv
charges.csv` writes the same rows to a CSV file instead, with charges at full precision, and
`--result-json` also writes the [result JSON](../guide/results.md). The CSV loads into analysis tools,
for example pandas:

```python
import pandas as pd

charges = pd.read_csv("charges.csv")
print(charges.groupby("element")["charge"].describe())
```

All records of the file are calculated with one method, parameter set, and execution mode. If you name
only the method, ChargeFW picks a parameter set that covers every record; if a record cannot be
covered, the calculation fails. The [Python page](../interfaces/python.md#one-model-for-all-molecules)
shows how to calculate records separately instead.

Download {download}`calculate_file.py`.

```{literalinclude} calculate_file.py
:language: python
```
