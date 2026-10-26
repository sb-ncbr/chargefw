# Bundled parameter sets

ChargeFW ships the parameter sets below. Within each method they are listed in the order automatic
selection tries them: ChargeFW uses the first set that covers every atom and bond of your molecules.
To get the same charges in every run, name the parameter set explicitly; see
[Methods and parameter sets](../guide/methods.md).

**Elements** lists the elements with atom parameters. **Atom types by** shows what besides the element
distinguishes atoms: the highest order of the atom's bonds, or the elements bonded to it. A set covers
an atom only if it has an entry for that element and type. The
[parameter-set format](parameter-format.md#classifiers) defines the types exactly.

On the command line, `chargefw parameters --method METHOD` lists the same sets; in Python, use
`chargefw.parameter_sets`.

```{include} ../_generated/parameter-sets.md
```
