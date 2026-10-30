# Command line

The `chargefw` command calculates charges for molecular files and writes the results to a directory.
The Python package installs it; [Installation](../getting-started/installation.md) also describes the
Docker image and building from source.

| Command | What it does |
| --- | --- |
| `chargefw calculate INPUT OUTPUT_DIRECTORY` | Calculates charges and writes the results |
| `chargefw inspect INPUT` | Reports what ChargeFW reads from a file |
| `chargefw applicability INPUT` | Reports which methods and parameter sets apply, without calculating |
| `chargefw methods [METHOD]` | Lists the methods, or describes one |
| `chargefw parameters [PARAMETER_SET]` | Lists the parameter sets, or describes one |

`chargefw COMMAND --help` lists the options of a command.

## Calculating charges

```console
$ chargefw calculate molecule.sdf output
Wrote output/molecule.chargefw.json
```

ChargeFW chooses the input format from the file extension: `.sdf` and `.mol`, `.mol2`, `.pdb`, `.cif`
and `.mmcif`, or `.json` for molecule JSON. It creates the output directory if needed and names the
results after the input file.

Without other options, ChargeFW chooses the method and parameter set itself, as described in
[Choosing a method](../guide/methods.md#choosing-a-method). Name them to get the same model every time:

```bash
chargefw calculate --method eem --parameter-set EEM_Cheminf_b3lyp_npa molecule.sdf output
```

ChargeFW calculates exactly what you request. If the method, parameter set, or execution mode cannot
handle the molecules, the calculation fails with the reasons; ChargeFW does not switch to another.

### Output files

Once ChargeFW has read the input and accepted the options, it writes `<input>.chargefw.json`, the
[result JSON](../guide/results.md) with all charges and the record of the calculation. It is written
also when the calculation fails or is cancelled, with the reasons. Errors before that point, such as a
malformed input file or an invalid option, are reported on standard error and leave no result file. On
request, ChargeFW also writes the molecules with their charges:

```bash
chargefw calculate --output-mol2 molecule.sdf output    # output/molecule.chargefw.mol2
chargefw calculate --output-mmcif structure.cif output  # output/structure.chargefw.cif
chargefw calculate --output-pqr structure.cif output    # output/structure.chargefw.pqr
```

MOL2 and mmCIF can be written for input in any format, but only for molecules with coordinates. PQR
is written only for a structure read from PDB or mmCIF, with one model; read a file with several
models using `--conformers first`. [Choosing an output
format](../guide/results.md#choosing-an-output-format) compares them.

ChargeFW writes the result JSON first. If writing MOL2, mmCIF, or PQR then fails, the JSON stays, still
with status `success`, and the command exits with status 6. A file that could not be written is not
created, and an existing file of that name from an earlier run stays as it was.

### Progress and cancelling

`--progress` shows the progress of the calculation. Press `Ctrl+C` once to request cancellation.
ChargeFW checks for it before each molecule, conformer, and fragment; when it stops there, it writes the
result JSON with status `cancelled` and no charges, and exits with status 5. A request that arrives
during the last calculation step comes too late, and the calculation finishes normally. Reading the
input and writing the results do not check for cancellation.

## Options

### Reading the input

| Option | Meaning |
| --- | --- |
| `--conformers first\|all` | Read only the first conformer or model, or all of them (default `all`) |
| `--structural-selection all\|polymers-and-ligands\|polymers` | PDB and mmCIF: which residues to keep; `polymers-and-ligands` leaves out water (default `all`) |
| `--structural-bonds none\|explicit\|templates\|hybrid` | PDB and mmCIF: where bonds come from (default `hybrid`, residue templates and bonds in the file) |

The [format reference](../reference/formats.md#pdb-and-mmcif-input) describes the selection and bond
options in detail. The structural options are an error with SDF, MOL, MOL2, or molecule JSON input.
A malformed record in the input stops the command before any calculation.

### Choosing the model

| Option | Meaning |
| --- | --- |
| `--method ID` | Use this method |
| `--parameter-set ID` | Use this parameter set; requires `--method` |
| `--method-option METHOD.OPTION=VALUE` | Set an option of a method, such as `peoe.iters=8`; repeatable |
| `--permissive-types` | Let a parameter set match atoms whose bond order has no exact entry, as described in [permissive matching](../reference/parameter-format.md#strict-and-permissive-matching) |
| `--fixed-ions ID` | Keep the ions of residue `ID` as fixed point charges, such as `MG`; repeatable; see [Fixed ions](../guide/fixed-ions.md) |

### Execution

| Option | Meaning |
| --- | --- |
| `--execution auto\|full\|cutoff\|cover` | Calculate the whole molecule or fragments (default `auto`); see [Large structures](../guide/large-structures.md) |
| `--radius ANGSTROM` | Fragment radius, at least 8 Å; required with `--execution cutoff` or `cover` (`auto` defaults to 12 Å) |
| `--cutoff-atom-threshold COUNT\|unlimited` | Size above which `auto` uses cutoff (default 20,000 atoms) |
| `--cover-atom-threshold COUNT\|unlimited` | Size above which `auto` uses cover (default 80,000 atoms) |
| `--threads COUNT` | Maximum number of threads (default `0`, all available cores) |

### Output and progress

| Option | Meaning |
| --- | --- |
| `--output-mol2` | Also write MOL2 |
| `--output-mmcif` | Also write mmCIF |
| `--output-pqr` | Also write PQR |
| `--progress` | Show the progress of the calculation |

`applicability` accepts the input, model, and execution options of `calculate`; `inspect` accepts the
input options.

## Checking before calculating

`chargefw inspect` reports what ChargeFW read, without loading any parameters:

```console
$ chargefw inspect 1l2y-h.cif
records: 1
record 0 (1L2Y) atoms=304 bonds=310 conformers=1 coordinates=true formal_charge=1 elements=H:150,C:98,N:27,O:29
```

[Preparing input](../guide/input.md#check-the-molecule-before-calculating) explains what to compare.

`chargefw applicability` lists every plan that could run, with its estimated memory and number of
operations, the reasons for every rejected method and parameter set, and the plan `calculate` would
choose. With `--fixed-ions` or execution options, it shows their effect without calculating.
[Large structures](../guide/large-structures.md#choosing-the-execution-yourself) shows an example.

## Methods and parameter sets

`chargefw methods` lists the methods; `chargefw methods METHOD` describes one, including its options:

```console
$ chargefw methods peoe
id: peoe
name: PEOE
full name: Partial Equalization of Orbital Electronegativity
publication: 10.1016/0040-4020(80)80168-2
notes: With the default zero-charge initialization, every connected component must be formally neutral. initial_charges=formal is a ChargeFW extension that preserves supplied formal charges.
requires coordinates: no
time complexity: O(n + m)
memory complexity: O(n)
supports cutoff: no
supports cover: no
supports fixed ions: no
options:
  iters (default=6) minimum>=1 maximum<=25
  initial_charges (default=zero) choices=zero,formal
```

`chargefw parameters` lists the parameter sets, `--method METHOD` only those of one method, and
`chargefw parameters ID` describes one. Add `--keys` to the last form to list every atom and bond type
the set has values for. The [parameter-set reference](../reference/parameter-sets.md) shows the same
sets with their publications.

## Exit statuses

| Status | Meaning |
| --- | --- |
| `0` | Success |
| `1` | Unexpected internal error |
| `2` | Invalid input or options, or the result JSON could not be written |
| `3` | No executable plan for the molecules and the requested settings |
| `4` | The calculation failed numerically, for example did not converge |
| `5` | The calculation was cancelled |
| `6` | Writing MOL2, mmCIF, or PQR failed after the result JSON was written |
