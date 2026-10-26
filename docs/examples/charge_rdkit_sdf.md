# Write charges to RDKit molecules

`charge_rdkit_sdf.py` reads an SDF with RDKit, calculates all records with one
method, parameter set, and execution mode, attaches the charges as the `ChargeFWPartialCharge` atom
property, and records the model as molecule properties. The output SDF stores the charges as an
`atom.dprop.ChargeFWPartialCharge` property list, which RDKit restores as atom properties when it reads
the file. Records that RDKit cannot read are reported and skipped. Every remaining record must convert,
so records need explicit hydrogens.

```bash
pip install "chargefw[rdkit]"
python charge_rdkit_sdf.py input.sdf charged.sdf --method eem
```

Because the converted records form one collection, ChargeFW selects a parameter set applicable to all of
them, and a record the selected model cannot handle makes the calculation fail rather than switching
models. Pass `--parameter-set` to choose the model explicitly. The example kekulizes aromatic rings so
ChargeFW receives explicit single and double bonds; see the
[RDKit integration](../interfaces/python.md#rdkit).

Download {download}`charge_rdkit_sdf.py`.

```{literalinclude} charge_rdkit_sdf.py
:language: python
```
