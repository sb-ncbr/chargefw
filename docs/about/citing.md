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
3.0.0 sqeqp 10.1186/s13321-021-00528-w
SQEqp_Schindler2021_CCD_gen 10.1186/s13321-021-00528-w
```

The publications are DOIs. Methods with built-in values, such as EQeq, MGC, and VEEM, have no
parameter set, and `plan.parameter_set` is then `None`. On the command line, `chargefw --version`
prints the version, and `chargefw methods ID` and `chargefw parameters ID` show the same DOIs, and the [parameter-set
reference](../reference/parameter-sets.md) lists them for every bundled set.

For ChargeFW itself, cite the publication of Atomic Charge Calculator III, whose calculation core
ChargeFW is ([Raček et al., 2026](https://doi.org/10.1093/nar/gkag379)), together with the ChargeFW
version you used. [CITATION.cff](https://github.com/sb-ncbr/chargefw/blob/main/CITATION.cff) has the
complete reference; on GitHub, **Cite this repository** formats it as APA or BibTeX.
