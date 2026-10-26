# Command-line interface

The `chargefw` executable imports molecular records, assesses available methods and parameter sets, and
writes source-mapped partial charges.

See [Installation](../getting-started/installation.md) for the Docker image and building from
source.

## Commands

```text
chargefw calculate [options] INPUT OUTPUT_DIRECTORY
chargefw inspect [input-options] INPUT
chargefw applicability [options] INPUT
chargefw methods [METHOD]
chargefw parameters [--method METHOD] [PARAMETER_SET]
```

Run `chargefw COMMAND --help` for the option syntax accepted by the installed version.

## Calculate charges

```bash
chargefw calculate molecule.sdf output
```

If no method or parameter set is specified, ChargeFW assesses the bundled catalog and chooses the first
plan in deterministic priority order. The output directory is created if necessary.

### Input options

Input format is selected from the file extension.

| Extensions | Reader |
| --- | --- |
| `.mol`, `.sdf` | Native MOL/SDF reader |
| `.mol2` | Native Tripos MOL2 reader |
| `.json` | ChargeFW molecule JSON 1.0 reader |
| `.pdb` | Gemmi-backed PDB reader |
| `.cif`, `.mmcif` | Gemmi-backed mmCIF reader |

The [molecular format reference](../reference/formats.md) describes the supported subsets, imported molecular data,
record mapping, warnings, and errors.

`--conformers first|all` selects the first conformer/model or all conformers/models. The default is
`all`.

PDB and mmCIF input additionally accepts:

- `--structural-selection all|polymers-and-ligands|polymers` (default `all`);
- `--structural-bonds none|explicit|templates|hybrid` (CLI default `hybrid`).

Structural options are rejected for other input formats. Selection, bonding, and alternate-location
semantics are defined in the [PDB and mmCIF format reference](../reference/formats.md#pdb-and-mmcif-input).

The CLI reads the imported collection before calculation and stops at the first malformed record.

### Method and execution options

These options are accepted by both `calculate` and `applicability`, except `--progress`, which applies
only to `calculate`. For `applicability`, `--threads` sets the thread limit used for plan estimates.

| Option | Meaning |
| --- | --- |
| `--method ID` | Restrict assessment to one method |
| `--parameter-set ID` | Parameter set for `--method` (required) |
| `--permissive-types` | Allow [permissive parameter classification](../reference/parameter-format.md#strict-and-permissive-matching) |
| `--method-option METHOD.OPTION=VALUE` | Override a method option; repeatable |
| `--fixed-ions ID` | Select an imported monatomic ion component with no graph bonds; repeatable |
| `--execution auto|full|cutoff|cover` | Select execution policy; default `auto` |
| `--radius ANGSTROM` | Radius for cutoff/cover; explicit reduced modes require at least 8 Å |
| `--cutoff-atom-threshold COUNT|unlimited` | Automatic full-to-cutoff threshold; default 20,000 |
| `--cover-atom-threshold COUNT|unlimited` | Automatic cutoff-to-cover threshold; default 80,000 |
| `--threads COUNT` | Maximum calculation threads; `0` delegates to oneTBB |
| `--progress` | Render calculation progress on standard error |

Use `chargefw methods` to list option IDs, defaults, allowed choices, and numeric bounds. Method options
are method-scoped even when a method is selected explicitly:

```bash
chargefw calculate --method peoe --method-option peoe.iters=8 molecule.sdf output
chargefw calculate --method sqeqp --parameter-set SQEqp_Schindler2021_CCD_gen \
    --fixed-ions MG --fixed-ions CA structure.cif output
chargefw applicability --fixed-ions MG --fixed-ions CA structure.cif
```

Named ions use the shared catalog and imported-component identity contract described in
[Fixed ions](../guide/fixed-ions.md). Successful calculations retain resolved
component IDs, source indices, and assigned charges in effective result provenance; see the
[result JSON reference](../reference/formats.md#chargefw-result-json-10).

Explicit method, parameter-set, or execution choices do not fall back to alternatives if they are
inapplicable. `applicability` reports the reasons without running a calculation.

Automatic and explicit execution share the policy, supported-method, conservation, and approximation
semantics documented in [Large structures](../guide/large-structures.md). Explicit choices
do not silently fall back to another mode.

## Output

The output basename is derived from the input filename:

```text
OUTPUT_DIRECTORY/<input-stem>.chargefw
```

Every normal calculation outcome writes result JSON:

```text
<basename>.json
```

On success, request fresh molecular representations explicitly:

```bash
chargefw calculate --output-mol2 molecule.sdf output
chargefw calculate --output-mmcif structure.cif output
chargefw calculate --output-mol2 --output-mmcif molecule.sdf output
```

The requested files use the same basename:

```text
<basename>.mol2
<basename>.cif
```

MOL2 and mmCIF are fresh generated representations for every supported input format and are not written
unless requested. The [molecular format reference](../reference/formats.md#charge-output) describes their intended use,
coordinate requirements, generated structures, charge fields, precision, and mapping checks.

Once import and request construction succeed, normal calculation outcomes write `<basename>.json`. On
success, ChargeFW writes JSON before attempting requested molecular exports. A molecular export failure
reports `Export error`, exits with status 6, and leaves JSON and any earlier completed export intact; the
JSON calculation status remains `success`. Files are written directly to their destinations, so a write
failure can leave that destination incomplete. Import, request-construction, primary JSON output, and other
filesystem failures are reported on standard error and can occur before a result document is published.

During calculation, press `Ctrl+C` once to request cooperative cancellation. ChargeFW stops at its next
[cancellation check point](../guide/results.md), writes `<basename>.json` with status
`cancelled` and no charge assignments, and exits with status 5. This does not interrupt import, request
construction, or output writing.

### JSON result

ChargeFW writes result JSON schema `1.0`; this is distinct from molecule input JSON schema `1.0`. See the
[result JSON format](../reference/formats.md#chargefw-result-json-10) for its records, diagnostics, assignments,
provenance, metrics, and numeric precision.

## Inspection and discovery

Inspect imported records without loading parameter sets:

```bash
chargefw inspect molecule.sdf
```

Report runnable plans, rejected alternatives, and the selected plan without calculating:

```bash
chargefw applicability --method eem --threads 4 molecule.sdf
```

Each runnable plan reports its estimated peak working memory and operation count for the thread limit;
the [resource estimates](../guide/large-structures.md#threads-and-resource-estimates) describes the estimate.

List method IDs and full names, inspect one method, or list bundled parameter sets:

```bash
chargefw methods
chargefw methods eem
chargefw parameters
chargefw parameters --method eem
chargefw parameters QEq_original
chargefw parameters SQEqp_Schindler2021_PUB_pept --keys
```

`chargefw methods METHOD` reports the selected method's identity, publication, human-readable notes,
priority, coordinate requirement, time and memory complexity, cutoff and cover support, fixed-point-source
support for `--fixed-ions`, and option schema.
Complexity uses Big-O notation with `n` for atoms and `m` for bonds, as defined in the
[Methods and parameter sets](../guide/methods.md).

`chargefw parameters` lists parameter-set IDs and names. Use `--method METHOD` to filter the summary by
the parameter set's declared method. Pass a parameter-set ID instead to report its method, name,
publication, notes, priority, and the number and classifiers of its atom and bond keys; the positional ID
and `--method` cannot be combined. Add `--keys` to a parameter-set ID to list every atom and bond key
in matching order, as `classifier=type` after the element symbol; an empty type is shown as `""`.
[Parameter-set JSON](../reference/parameter-format.md#classifiers) defines classifiers and matching precedence.

## Exit statuses

| Status | Meaning |
| --- | --- |
| `0` | Calculation or reporting command completed successfully |
| `1` | Unexpected internal failure |
| `2` | Invalid input or request, or primary output publication failure |
| `3` | No executable plan |
| `4` | Numerical calculation failure |
| `5` | Calculation cancelled |
| `6` | A requested molecular export failed after JSON publication |
