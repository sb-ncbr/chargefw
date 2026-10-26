# Inspect imported molecules

`inspect_molecules.py` reads any supported molecular format and reports the source
identity, dimensions, conformers, and total formal charge of every imported record. Structural selection,
bond, and conformer policies are explicit options. With `--classifier plain`, `hbo`, or `bonded`, it also
counts each record's atoms by element and [parameter classifier](../reference/parameter-format.md#classifiers) type,
showing how finely that classifier separates the atoms, for example before grouping calculated charges.

```bash
python inspect_molecules.py input.sdf --format sdf --classifier bonded
```

Download {download}`inspect_molecules.py`.

```{literalinclude} inspect_molecules.py
:language: python
```
