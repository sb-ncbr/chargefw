# Inspect imported molecules

This script shows what ChargeFW reads from a file: for every record, its identity, atom and bond counts,
conformers, and total formal charge. With `--classifier`, it also counts the atoms by element and by
the [atom type](../reference/parameter-format.md#classifiers) that parameter sets use, which shows which
atom types a parameter set needs to cover the molecule:

```console
$ python inspect_molecules.py paracetamol.sdf --format sdf --classifier hbo
Collection: paracetamol.sdf (1 molecule(s))
[0] record=0 id='paracetamol' name='paracetamol'
    atoms=20 bonds=20 conformers=1 formal_charge=0
    conformer_names=('input',)
    C   hbo=1: 1
    C   hbo=2: 7
    H   hbo=1: 9
    N   hbo=1: 1
    O   hbo=1: 1
    O   hbo=2: 1
```

`hbo` distinguishes atoms by their highest bond order; `plain` only by element and `bonded` by their
bonded neighbors. For PDB and mmCIF, the selection, bond, and conformer options of
[reading molecules](../interfaces/python.md#reading-molecules) are available as command-line options.

Download {download}`inspect_molecules.py`.

```{literalinclude} inspect_molecules.py
:language: python
```
