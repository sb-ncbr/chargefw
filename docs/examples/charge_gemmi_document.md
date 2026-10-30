# Charge a Gemmi document

This script reads an mmCIF file with Gemmi, calculates charges, and adds them to the document. All
other data of the file is kept. The output is mmCIF with the charges in the SB NCBR partial atomic
charges categories. For the Trp-cage miniprotein [1l2y-h.cif](../guide/1l2y-h.cif):

```console
$ pip install --pre "chargefw[gemmi]"
$ python charge_gemmi_document.py 1l2y-h.cif charged.cif
```

`charged.cif` then contains:

```text
_sb_ncbr_partial_atomic_charges.charge
1 1 -0.8727
1 2 -0.0952
1 3 0.5663
...
```

The script uses EEM with a parameter set ChargeFW chooses. `1l2y-h.cif` is prepared, with hydrogens
and formal charges; your own structures need the same, as described in
[Preparing input](../guide/input.md). The
[Gemmi section](../interfaces/python.md#gemmi) of the Python page describes the checks ChargeFW makes
before it changes the document.

Download {download}`charge_gemmi_document.py`.

```{literalinclude} charge_gemmi_document.py
:language: python
```
