# Analyze conformer-dependent charges

`analyze_rdkit_conformer_charges.py` generates an unoptimized,
deterministic RDKit ETKDG ensemble, calculates EEM charges, and ranks atoms by charge standard deviation
over the equally weighted generated conformers.

```bash
pip install "chargefw[rdkit]"
python analyze_rdkit_conformer_charges.py "CCO" --conformers 20
```

This is a sensitivity analysis over one generated coordinate set, not a thermodynamic ensemble or a
claim that the ranked atoms are chemically significant hotspots. The useful analysis function also
accepts an externally prepared RDKit conformer ensemble.

Download {download}`analyze_rdkit_conformer_charges.py`.

```{literalinclude} analyze_rdkit_conformer_charges.py
:language: python
```
