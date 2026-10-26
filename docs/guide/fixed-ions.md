# Fixed ions

Many structures contain monatomic ions such as Mg²⁺, Zn²⁺, Na⁺, or Cl⁻. Most parameter sets have no
values for them, so a method that has to include them in its calculation cannot be used. Removing the
ions from the structure avoids that, but their charge then has no effect on the rest of the molecule.

With fixed ions, ChargeFW keeps each selected ion at its formal charge, as a point charge, and calculates
the charges of all other atoms in the electric field of the ions. EEM, SQE, SQE+q0, and SQE+qp support
fixed ions; the **Fixed ions** column of the [method table](methods.md#available-methods) shows them.

## Selecting ions

Select ions by their residue name in the structure, such as `MG` or `CA`. On the command line, repeat
`--fixed-ions`:

```console
$ chargefw calculate --fixed-ions MG --fixed-ions CA structure.cif output
```

In Python, pass the names to `fixed_ions`, or use `chargefw.COMMON_IONS` for Na⁺, K⁺, Mg²⁺, Ca²⁺,
Cl⁻, Zn²⁺, Fe²⁺, and Fe³⁺, or `chargefw.ALL_IONS` for the whole [catalog](#ion-catalog):

```python
result = chargefw.calculate(molecules, fixed_ions=["MG", "CA"])
result = chargefw.calculate(molecules, fixed_ions=chargefw.COMMON_IONS)
```

A selected residue must consist of a single atom of the ion's element, without bonds. Names of ions that
are not in the structure are ignored; a name that is not in the catalog is an error.

Fixed ions work with PDB and mmCIF input, where atoms belong to named residues. Ions that are not
selected stay part of the calculated structure and need parameters like any other atom.

## What the result contains

In {download}`water-ions.cif`, a water molecule sits between a Mg²⁺ and a Ca²⁺ ion:

```console
$ chargefw calculate --method sqeqp --parameter-set SQEqp_Schindler2021_CCD_gen \
    --fixed-ions MG --fixed-ions CA water-ions.cif output
Wrote output/water-ions.chargefw.json
```

The charges in atom order are −0.877 for the oxygen, 0.422 and 0.455 for the hydrogens, and 2.000 for
each ion; their sum is 4.0. The water is polarized by the ions, and each ion has the charge of its
catalog entry. The charge comes from the catalog even when the file gives the ion a different formal
charge; the two ions here have formal charge 0 in the file. The result records which ions were fixed:

```json
"fixed_ions": [
  {"component_id": "MG", "charge": 2.0, "instances": [{"molecule_index": 0, "atom_index": 3}]},
  {"component_id": "CA", "charge": 2.0, "instances": [{"molecule_index": 0, "atom_index": 4}]}
]
```

In Python, `result.plan.fixed_ions` holds the same information.

The total charge of the calculated atoms is the sum of their formal charges; the fixed ions add their
own charges to it. The total of the whole result can therefore differ from the formal charges in the
file when the file gives an ion another charge than its catalog entry.

## What the approximation does

A fixed ion is a point charge that does not change. Its electric field can redistribute the charges of
the other atoms within the charge constraints of the method. The ion itself is not polarized, and no
charge moves to or from it. This approximation has not been validated for quantitative predictions of
metal coordination.

In EEM, the potential of the fixed ions enters the electronegativity of every calculated atom $i$:

$$
\chi_i = A_i + B_i q_i + \kappa \sum_{j \neq i} \frac{q_j}{R_{i,j}} + \kappa \sum_F \frac{q_F}{R_{i,F}}
$$

where $q_F$ is the charge of fixed ion $F$ and $\kappa$ is the same parameter as for the other atoms.

In SQE, SQE+q0, and SQE+qp, each ion adds the potential

$$
\varphi_{i,F} = q_F \, \frac{\operatorname{erf}\left(\dfrac{R_{i,F}}{\sqrt{2}\,|w_i|}\right)}{R_{i,F}}
$$

to the electronegativity of atom $i$, where $w_i$ is the Gaussian width of atom $i$; with zero width, the
potential is $q_F / R_{i,F}$. The [method descriptions](methods.md#method-descriptions) give the full
equations.

With [cutoff and cover](large-structures.md), every fragment contains all fixed ions, including those
beyond the radius, and the ions never become fragment centers.

## Ion catalog

ChargeFW recognizes 80 monatomic ions from the Chemical Component Dictionary, each with its charge:

| Charge | Residue names |
| --- | --- |
| −1 | `BR`, `CL`, `F`, `IOD` |
| +1 | `AG`, `AU`, `CS`, `CU1`, `D8U`, `K`, `LI`, `NA`, `RB`, `RH`, `TL` |
| +2 | `0BE`, `BA`, `CA`, `CD`, `CO`, `CU`, `EU`, `FE2`, `HG`, `MG`, `MN`, `NI`, `PB`, `PD`, `PT`, `RHF`, `SR`, `Y1`, `YB2`, `ZN` |
| +3 | `3CO`, `3NI`, `AL`, `AM`, `AU3`, `BS3`, `CE`, `CF`, `CR`, `CU3`, `DY`, `ER3`, `EU3`, `FE`, `GA`, `GD3`, `HO3`, `IN`, `IR3`, `LA`, `LU`, `MN3`, `ND`, `OS`, `PR`, `RH3`, `RU`, `SB`, `SM`, `TB`, `V`, `YB`, `YT3`, `ZCM`, `ZTM` |
| +4 | `4MO`, `4PU`, `4TI`, `IR`, `OS4`, `PT4`, `TH`, `ZR` |
| +6 | `6MO`, `W` |

`CA` is calcium and `FE` iron(III), as in the Chemical Component Dictionary. The
[fixed-ion example](../examples/calculate_with_fixed_ions.md) reports the charges around every fixed ion
of a structure.
