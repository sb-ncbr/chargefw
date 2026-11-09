# Fitting parameters

`chargefw.parameterisation` fits the parameter values of a method to reference charges you already
have, for example charges from quantum chemistry for a set of training molecules. The result is a
parameter set you can calculate with and save like any other. Fitting needs SciPy 1.18 or newer;
install it with the `parameterisation` extra:

```console
$ pip install "chargefw[parameterisation]"
```

This fits SQE+qp from scratch to reference charges stored as one text file per molecule, with one
charge per line in atom order:

```python
import numpy as np
import chargefw
from chargefw.parameterisation import Evaluator, fit, make_parameter_template

molecules = chargefw.io.read("training.sdf", format="sdf")
references = [np.loadtxt(f"charges/{index}.txt") for index in range(len(molecules))]

template = make_parameter_template(molecules, "sqeqp", atom_types="hbo", bond_types="bo")
plan = chargefw.assess(
    molecules, method="sqeqp", parameter_set=template, execution="full"
).default_plan
evaluator = Evaluator(plan)

# Electronegativities start from the Pauling values of the elements in your molecules.
pauling = {"H": 2.20, "C": 2.55, "N": 3.04, "O": 3.44, "S": 2.58}
start = {
    (slot.key, slot.name): pauling[slot.key.element]
    for slot in evaluator.layout
    if slot.name == "electronegativity"
}

result = fit(evaluator, references, start=start)
print(result.message, result.rmsd, result.diagnostics.rank, len(result.free))

fitted = result.parameter_set(id="my-sqeqp", name="SQE+qp for my molecules")
chargefw.write_parameter_set(fitted, "my-sqeqp.json")
```

The rest of this page explains each step and the choices you can make along the way. The examples
continue from this one.

## Making a template

A parameter set has one entry per atom type and bond type, and atoms of the same type share their
parameter values. `make_parameter_template(molecules, method, atom_types=..., bond_types="plain")`
finds the types that occur in your molecules and makes a set with an entry for each of them, all
values zero. The types come from the classifiers of the [parameter-set
format](../reference/parameter-format.md):

| `atom_types` | Atoms of one type share | For example, oxygen in |
| --- | --- | --- |
| `"plain"` | the element | every molecule is one type |
| `"hbo"` | the element and the highest order of its bonds | an alcohol (`O/1`) and a carbonyl (`O/2`) are two types |
| `"bonded"` | the element and the elements of its bonded atoms | an alcohol (`O/CH`) and an ether (`O/CC`) are two types too |

`bond_types` is `"plain"` (one type for each pair of atom types) or `"bo"` (also split by bond order).

Finer types can describe more chemistry, but each type has its own values, and those need enough
atoms of that type in your molecules to be determined. With `"bonded"`, a type that occurs in only one
or two molecules is common; its values then fit those few atoms rather than the chemistry. Start with
`"hbo"` and use finer types when you have many training molecules; [reading the
result](#reading-the-result) shows which types your data cannot pin down.

You do not have to start from zero. To improve a parameter set you fitted before, for example with
more molecules, assess your molecules with it and fit from its values:

```python
previous = chargefw.read_parameter_set("my-sqeqp.json")
plan = chargefw.assess(
    molecules, method=previous.method, parameter_set=previous, execution="full"
).default_plan
evaluator = Evaluator(plan)
result = fit(evaluator, references)  # starts from the set's own values
```

Every atom and bond of your molecules must then match an entry of that set. A bundled or published
set works the same way, but its values must satisfy the profile's limits: some published SQE sets
have a hardness below the stability limit described below, and `fit()` then names those values and
stops before fitting.

## Calculating with candidate values

`Evaluator(plan)` calculates the charges of your molecules for any parameter values. It reads and
prepares the molecules once, so a fit can try thousands of values quickly. The plan must use full
execution (`execution="full"` in `assess()`), so that the charges are exact and not affected by
[cutoff or cover execution](../guide/large-structures.md).

The values are a vector in a fixed order, `evaluator.layout`. Each entry, a `Slot`, names one value:
its `scope` (`"common"`, `"atom"`, or `"bond"`), the `key` of its atom or bond type (`None` for a
common parameter), and the parameter `name`. Print it to see what you are fitting:

```python
for index, (slot, count) in enumerate(zip(evaluator.layout, evaluator.matches)):
    print(index, slot.scope, slot.key, slot.name, count)
```

`evaluator.matches` counts the atoms or bonds of each type in your molecules: a type with only a few
is poorly determined by your data. `evaluator.charges(values)` returns the charges for a complete
vector of values, all assignments concatenated. A template's own values are all zero, which most
methods cannot calculate; after a fit, `result.start` holds the complete start the fit used and
`result.values` the fitted one, so you can compare the charges of both:

```python
start_charges = evaluator.charges(result.start)
fitted_charges = evaluator.charges(result.values)
```

## Fitting

`fit(evaluator, references, ...)` changes the values until the charges are as close to your
references as it can make them, measured by the sum of squared differences. Give `references` as
one array per assignment, as in the first example, or as one array for all charges. A method that
uses coordinates calculates every conformer, so it then needs one array per conformer.

The other arguments let you decide what to fit and how much each charge counts. You can leave them
all out: for EEM and the SQE methods, the method's *profile* then provides sensible defaults.

### What the profile does

| Method | Fixed | Limits | Start values |
| --- | --- | --- | --- |
| EEM | `kappa`, and `A` of the first type | `B ≥ 0.001` | `kappa = 1`, `B = 2` |
| SQE, SQE+q0 | `electronegativity` of the first type, and the `hardness` of terminal types | `width ≥ 0.05`, `kappa ≥ 0.01`, a stable hardness | `hardness = 2`, `width = 1`, `kappa = 1` |
| SQE+qp | as SQE, and `q0` of the first type | as SQE | as SQE, and `q0 = 0` |

Some changes of the values do not change any charge. Adding the same amount to every
electronegativity is one: only differences between electronegativities move charge. Such a direction
cannot be fitted, so the profile fixes one value to remove it. The others are fitted relative to that
one.

In the SQE methods, a *terminal* type, one whose atoms all have exactly one bond in your molecules
(hydrogen, usually halogens, often a carbonyl oxygen), has another: such an atom exchanges charge only
through its bond, so its hardness and the `kappa` of its bonds act only as a sum. When those bonds
have parameter entries of their own, as in a template, the profile holds the hardness of the terminal
type at its stability limit and raises the `kappa` of its bonds by the difference, which leaves every
charge unchanged, also of your start. This costs no accuracy, and it makes the fitted values
reproducible. The width of a terminal type is still fitted.

The limits keep the values where the method works. In the SQE methods, too small a hardness lets
charge flow without limit; the profile fits each hardness so that √π·hardness·width ≥ 1, which keeps
the calculation stable for every molecule and geometry, also for molecules you calculate later.

Electronegativities have no default start. Give them in `start`, with different values for different
elements, such as the Pauling electronegativities in the first example. With equal
electronegativities, neutral molecules get zero charges, and no other parameter changes that, so the
fit could not start.

To set everything yourself, pass `profile=None`.

### Choosing what to fit

`free` lists the slots to fit, as indices into `evaluator.layout`; the other slots keep their start
values. Use it when you trust some values and want to fit only the rest. For example, to add sulfur
to a parameter set that already describes your other elements, fit only the sulfur types and their
bonds:

```python
free = [
    index
    for index, slot in enumerate(evaluator.layout)
    if (slot.scope == "atom" and slot.key.element == "S")
    or (slot.scope == "bond" and any(atom.element == "S" for atom in slot.key.atoms))
]
result = fit(evaluator, references, free=free)
```

Fixing values also helps when the data cannot determine them, as [reading the
result](#reading-the-result) explains.

In the SQE methods, the profile fits each hardness together with the width of the same type, to keep
the calculation stable. If you fit a width, fit its hardness too; `fit()` names the slots otherwise.

### Weighting charges

By default every charge counts equally. `weights` changes that: give one weight per charge, in the
same shape as `references`. A charge with weight 2 counts as much as four with weight 1, because the
squared differences are weighted by the squares of the weights. Weights help when equal counting
does not fit your purpose:

- **Elements in very different numbers.** Hydrogen is often half of all atoms, so a fit with equal
  weights cares mostly about hydrogen charges. To give every element the same total influence,
  weight each charge by one over the square root of the number of charges of its element.
- **References of different quality.** Give charges you trust less a smaller weight.
- **Molecules to leave out.** Weight 0 removes a charge from the fit, but it is still calculated and
  reported. This keeps the template and the plan unchanged, for example when you compare fits on
  different subsets.

```python
elements = [
    np.array([chargefw.periodic_table[z].symbol for z in molecule.atomic_numbers])
    for molecule in molecules
]
counts = {}
for symbols in elements:
    for symbol in symbols:
        counts[symbol] = counts.get(symbol, 0) + 1
weights = [np.array([counts[symbol] ** -0.5 for symbol in symbols]) for symbols in elements]
result = fit(evaluator, references, start=start, weights=weights)
```

`result.rmsd` and the errors in the [diagnostics](#reading-the-result) are reported without weights,
so you can compare fits with different weights.

### Limits

`bounds` gives a lower and an upper value per slot, as two arrays in layout order; use `-np.inf` and
`np.inf` for no limit. Use it to keep a value in a range you consider meaningful, for example a
positive `kappa`. Bounds you give replace the profile's, with one exception: in the SQE methods, the
stability limit of each hardness stays, and a hardness takes no bounds of its own. To set every limit
yourself, pass `profile=None`; the fit then no longer guarantees stable values.

A value that ends at its limit is reported in the [diagnostics](#reading-the-result): the fit wanted
to go further, so check whether the limit or the type is the problem.

### Trying a fit quickly

`max_evaluations` limits how many trial values the optimizer may try. A small number, such as 20,
gives a quick first result, which shows whether the start, the types, and the weights make sense
before you spend time on a full fit. Calculations for checking the start, for the derivatives, and
for the diagnostics come on top of that number.

### Checks and failures

Before it starts, `fit()` checks that the start lies within the limits (a value on a limit, as an
earlier fit can leave it, is fine), that the charges can be calculated, and that every value you fit
changes at least one charge. Otherwise it raises `ValueError` naming the slots, or
`chargefw.parameterisation.EvaluationError` if the charges cannot be calculated.

During the fit, values for which the calculation fails, for example because the system becomes
unstable, are rejected, and the fit tries a smaller step. `result.failed_evaluations` counts them.

For EEM, SQE, SQE+q0, and SQE+qp, `fit()` uses exact derivatives of the charges: each step costs
about as much as a few calculations of the charges, however many values you fit. Other methods use
finite differences, which cost one calculation per fitted value and step, so they are slow with many
types.

## Reading the result

`fit()` returns a `FitResult`. Start with how it ended and how close the charges are:

```python
print(result.success, result.message, result.evaluations)
print(f"RMSD {result.rmsd:.4f} e, R² {result.diagnostics.r2:.4f}")
for symbol, element in result.diagnostics.by_element.items():
    print(f"{symbol}: {element.count} charges, RMSD {element.rmsd:.4f} e")
```

`result.rmsd` and `r2` compare the charges at the fitted values with your references, without weights.
`by_element` shows whether one element is fitted much worse than the others, which often means it
needs finer types, more molecules, or a larger weight.

Then check whether the values can be trusted. A good fit of the charges does not guarantee
meaningful values:

- **`rank` below the number of fitted values.** Some combination of values changes the charges so
  little that your data cannot determine it. For EEM and the SQE methods, the fit does not try to
  improve the charges through combinations like this: that would only work by pushing the values to
  extreme numbers. The values involved can still change during the fit, sometimes a lot, but the data
  does not decide them: another start can give quite different values with almost the same charges,
  so do not read their exact numbers as results. `weakest` names them, and `diagnostics.matches` shows how many atoms or bonds their types
  have. A typical cause is a type that occurs once, or, in the SQE methods, a terminal type whose
  hardness you fit yourself with `free`. To determine them, fix one of the values involved with `free`, use coarser
  types, or add training molecules in which those types occur in different surroundings.
- **`at_bounds`.** These values ended at a limit. The fit wanted to go beyond it, so the limit is
  shaping the result.
- **`failed_messages`.** Messages of calculations that failed during the fit, at most 20. A few are
  normal; many mean the fit spent its time near values the method cannot use.

`singular_values` show how strongly the data determines each independent combination of the fitted
values, strongest first; `rank` counts those that are not negligible.

Finally, judge the fit on molecules it did not see. Calculate a separate test set with the fitted
values and compare the charges with their references, per element and per molecule, and look at the
molecules with the largest errors. A fit that is much better on its training molecules than on the
test set has too many types for its data.

## Using the fitted values

`result.parameter_set(id=..., name="", publication="", notes="")` turns the fitted values into a
parameter set for the evaluator's method. The ID must differ from the IDs of the bundled sets. Use it
like any parameter set, for example on your test molecules:

```python
test = chargefw.io.read("test.sdf", format="sdf")
calculated = chargefw.calculate(test, method=fitted.method, parameter_set=fitted)
```

`chargefw.write_parameter_set(fitted, path)` saves it in the [parameter-set
format](../reference/parameter-format.md), and `chargefw.read_parameter_set(path)` reads it back.
