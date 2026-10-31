# Citing ChargeFW

When you publish charges calculated with ChargeFW, cite the method and the parameter set you used,
and ChargeFW itself. Report the method ID, the parameter-set ID, and the ChargeFW version, so that
others can repeat the calculation.

The [result of a calculation](../getting-started/first-calculation.md#python) names the model that
was used and the publications behind it:

```python
plan = result.plan
print(chargefw.__version__, plan.method.id, plan.method.publication)
print(plan.parameter_set.id, plan.parameter_set.publication)
```

```text
3.0.0rc2 sqeqp 10.1186/s13321-021-00528-w
SQEqp_Schindler2021_CCD_gen 10.1186/s13321-021-00528-w
```

The publications are DOIs. Methods with built-in values, such as EQeq, MGC, and VEEM, have no
parameter set, and `plan.parameter_set` is then `None`. On the command line, `chargefw methods ID`
and `chargefw parameters ID` show the same DOIs, and the [parameter-set
reference](../reference/parameter-sets.md) lists them for every bundled set.

For ChargeFW itself, cite the software with the metadata in
[CITATION.cff](https://github.com/sb-ncbr/chargefw/blob/main/CITATION.cff), together with the
version you used. On GitHub, **Cite this repository** formats it as APA or BibTeX.
