# Compare bundled parameter sets

This script calculates one molecule with every parameter set of a method that applies to it, and
compares the charges with those of the first set. For
[paracetamol.sdf](../getting-started/paracetamol.sdf):

```console
$ python compare_parameter_sets.py paracetamol.sdf --format sdf --method eem
Reference: EEM_NEEMP_ccd2016_npa
Parameter set                              total charge    RMS difference    max difference
EEM_NEEMP_ccd2016_npa                         -0.000000          0.000000          0.000000
EEM_NEEMP_ccd2016_npa2                         0.000000          0.044750          0.133496
EEM_NEEMP_ccd2016_mpa                          0.000000          0.081309          0.190372
EEM_NEEMP_ccd2016_mpa2                        -0.000000          0.067589          0.220710
EEM_Cheminf_b3lyp_aim                         -0.000000          0.228376          0.503868
...
```

The script [assesses](../interfaces/python.md#assessing-before-calculating) the molecule once and then
calculates each plan, so the molecule is prepared only once. It also lists the parameter sets that do not
apply, with their reasons.

The differences show how much the charges depend on the parameter set. Parameter sets were fitted to
different reference charges and molecules, so the comparison does not tell which one is more
accurate.

Download {download}`compare_parameter_sets.py`.

```{literalinclude} compare_parameter_sets.py
:language: python
```
