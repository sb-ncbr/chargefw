# Charge a Gemmi document

`charge_gemmi_document.py` attaches charges to the unchanged live Gemmi mmCIF
document used for import. PDB input is first represented as a caller-owned Gemmi mmCIF document. See the
[Gemmi integration](../interfaces/python.md#gemmi) for validation and mutation behavior.

```bash
pip install "chargefw[gemmi]"
python charge_gemmi_document.py input.cif charged.cif --format mmcif
```

The example uses EEM by default and lets ChargeFW select an applicable bundled parameter set.

Download {download}`charge_gemmi_document.py`.

```{literalinclude} charge_gemmi_document.py
:language: python
```
