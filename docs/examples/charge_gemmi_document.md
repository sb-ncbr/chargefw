# Charge a Gemmi document

This script reads an mmCIF or PDB file with Gemmi, calculates charges, and adds them to the document.
For mmCIF input, all other data of the file is kept. A PDB file is first converted to a new mmCIF
document, so only what Gemmi takes over from the PDB file is kept. The output is mmCIF with the charges
in the SB NCBR partial atomic charges categories:

```console
$ pip install --pre "chargefw[gemmi]"
$ python charge_gemmi_document.py structure.cif charged.cif --format mmcif
```

```text
_sb_ncbr_partial_atomic_charges.charge
1 1 -0.8661
1 2 -0.1166
1 3 0.5803
...
```

The script uses EEM with a parameter set ChargeFW chooses. The structure must be prepared, with
hydrogens and formal charges, as described in [Preparing input](../guide/input.md). The
[Gemmi section](../interfaces/python.md#gemmi) of the Python page describes the checks ChargeFW makes
before it changes the document.

Download {download}`charge_gemmi_document.py`.

```{literalinclude} charge_gemmi_document.py
:language: python
```
