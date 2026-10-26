# Calculate a molecular file

`calculate_file.py` is the shortest complete file-to-charges workflow. It reads any
supported molecular format, calculates one explicitly selected method, and prints one row per atom in
source order with its molecule and conformer indices, atom index, atom ID, atom name, element, and
charge. A `-` conformer identifies charges from a geometry-independent method. `--csv` writes the same
rows to a CSV file instead, with full-precision charges and an empty conformer field for
geometry-independent methods. `--result-json` also writes complete result JSON with requested and
effective provenance.

```bash
python calculate_file.py molecules.sdf --format sdf --method eem --csv charges.csv
```

The CSV loads directly into analysis tools; for example, with pandas installed:

```python
import pandas as pd

charges = pd.read_csv("charges.csv")
print(charges.groupby("element")["charge"].describe())
```

One method, parameter set, and execution mode apply to the whole file, so one inapplicable record makes
the request inapplicable; the example does not silently choose a separate model per record. Its EEM
invocation lets ChargeFW select a bundled parameter set applicable to every record. The
[Python collection guide](../interfaces/python.md#reading-molecular-files) explains when to calculate a collection or
its records independently, and the [format reference](../reference/formats.md) defines import and result JSON
behavior.

Download {download}`calculate_file.py`.

```{literalinclude} calculate_file.py
:language: python
```
