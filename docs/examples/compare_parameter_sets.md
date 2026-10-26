# Compare bundled parameter sets

This script calculates one molecule with every parameter set of a method that applies to it, and
compares the charges with those of the first set:

```console
$ python compare_parameter_sets.py paracetamol.sdf --format sdf --method eem
Reference: EEM_Ionescu2013_mpa_gas
Parameter set                              total charge    RMS difference    max difference
EEM_Ionescu2013_mpa_gas                        0.000000          0.000000          0.000000
EEM_Ionescu2013_mpa_pcm                        0.000000          0.030521          0.069620
EEM_Ionescu2013_npa_gas                       -0.000000          0.056624          0.129664
EEM_Ionescu2013_npa_pcm                       -0.000000          0.064932          0.158402
EEM_Baek1991                                   0.000000          0.237594          0.495312
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
