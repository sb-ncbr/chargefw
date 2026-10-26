# First calculation

This page calculates charges for paracetamol, first on the command line and then in Python. Download
{download}`paracetamol.sdf`, a structure with all hydrogens and 3D coordinates, to follow along.

## Command line

Give `chargefw calculate` the input file and a directory for the results:

```console
$ chargefw calculate paracetamol.sdf output
Wrote output/paracetamol.chargefw.json
```

With the [Docker image](installation.md#command-line-in-docker), run the same command through
`docker run` with your directory mounted.

You did not choose a method, so ChargeFW chose one: it tried the methods in their
[selection order](../guide/methods.md#available-methods) and used the first one with a parameter set
covering every atom of paracetamol. The result file records what was used:

```json
"calculation_provenance": {
  "effective": {
    "method": {"id": "sqeqp"},
    "parameter_set": {"id": "SQEqp_Schindler2021_CCD_gen"},
    "execution": {"mode": "full", "radius_angstrom": null}
  }
}
```

The charges are in `results[0].assignments[0].charges`, one value in elementary charges per atom, in
the order of the atoms in the input file: −0.637 for the first carbon, 0.643 for the carbonyl carbon,
and so on.

To choose the model yourself, name the method and the parameter set. `--output-mol2` also writes the
molecule with its charges in MOL2, which molecular viewers and docking programs read:

```console
$ chargefw calculate --method eem --parameter-set EEM_Cheminf_b3lyp_npa --output-mol2 paracetamol.sdf output
Wrote output/paracetamol.chargefw.json
Wrote output/paracetamol.chargefw.mol2
```

The charges are in the last column of the MOL2 atom records, rounded to four decimal places:

```text
@<TRIPOS>ATOM
1 C -3.6645 -0.5774 -0.0163 C 1 UNL -0.6159
2 C -2.2115 -0.8955 -0.2688 C 1 UNL 0.5291
3 O -1.8819 -1.9716 -0.7542 O 1 UNL -0.6475
```

## Python

The same calculation in Python:

```python
import chargefw

molecules = chargefw.io.read("paracetamol.sdf", format="sdf")
result = chargefw.calculate(molecules)

print(result.plan.method.id, result.plan.parameter_set.id)
print(result.assignments[0].values.round(3))
```

```text
sqeqp SQEqp_Schindler2021_CCD_gen
[-0.637  0.643 -0.607 -0.576  0.154 -0.168 -0.268  0.289 -0.696 -0.231
 -0.238  0.22   0.233  0.21   0.38   0.215  0.206  0.459  0.21   0.202]
```

`values` is a NumPy array in input atom order. Pass `method="eem"` and
`parameter_set="EEM_Cheminf_b3lyp_npa"` to `calculate()` to choose the model, as on the command line.

## Next steps

- [Preparing input](../guide/input.md): what your own molecules need before a calculation.
- [Methods and parameter sets](../guide/methods.md): which model to choose.
- [Results](../guide/results.md): conformers, atom order, and what the result records.
- [Examples](../examples/index.md): complete scripts for files, RDKit, Gemmi, and analysis.
