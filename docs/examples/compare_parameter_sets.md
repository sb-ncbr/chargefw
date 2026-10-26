# Compare bundled parameter sets

`compare_parameter_sets.py` performs one explicit assessment for a method
and molecular geometry, reports rejected parameter sets, and executes each prepared plan to compare its
source-aligned charge vector against one reference. This advanced workflow avoids repeating molecule
preparation and parameter classification for every executable parameter set.

```bash
python compare_parameter_sets.py molecule.sdf --format sdf --method eem
```

Differences demonstrate parameter sensitivity only. Parameter sets may target different molecular
domains or reference charge schemes; the example does not rank their scientific accuracy.

Download {download}`compare_parameter_sets.py`.

```{literalinclude} compare_parameter_sets.py
:language: python
```
