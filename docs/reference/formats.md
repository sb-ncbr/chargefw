# Molecular input and charge output formats

This document describes the format behavior shared by ChargeFW's native adapters, command-line
application, and language bindings. API names and command-line file-selection policy are documented in
[Native C++ library](../interfaces/cpp.md), [Python package](../interfaces/python.md), and
[Command-line interface](../interfaces/cli.md), respectively.

ChargeFW imports files into a toolkit-neutral molecule containing source-ordered atoms, bonds, formal
charges, and zero or more Cartesian conformers. Readers apply only the format-specific transformations
documented below; broader chemistry-preparation limits are defined in the
[Preparing input](../guide/input.md). Atom and molecule order is retained in charge
assignments and output mappings.

## Format overview

| Format | Accepted input | Available output | Input record model |
| --- | --- | --- | --- |
| MOL V2000/V3000 | Yes | No | One molecule, one conformer |
| SDF V2000/V3000 | Yes | No | One molecule per record, one conformer each |
| Tripos MOL2 | Yes | Yes, freshly generated | One molecule per `MOLECULE` record, one conformer each |
| ChargeFW molecule JSON 1.0 | Yes | No | `molecules` array, zero or more conformers each |
| PDB | Yes, through Gemmi | No | One molecule, models become conformers |
| mmCIF | Yes, through Gemmi | Yes, freshly generated | One molecule per coordinate-bearing block |
| ChargeFW result JSON 1.0 | No | Yes | Not applicable |

Input and output capabilities are independent. A result calculated from any supported molecular input can
be written as result JSON or, when its data is representable, as generated MOL2 or generated mmCIF. ChargeFW
does not write MOL, SDF, molecule JSON, or PDB.

Readers report malformed or unsupported molecular data as errors. Native MOL, SDF, and MOL2 input accepts
both LF and CRLF line endings. A successful imported record also carries its source name, zero-based record
index, format-derived record ID, and non-fatal diagnostics. MOL, SDF, MOL2, and molecule JSON records own
their calculation-atom source positions, explicit source atom IDs where available, per-conformer site
mapping, import policy, and whether source connectivity was absent, explicitly empty, or present. Streaming
SDF and MOL2 readers consume one record at a time; molecule JSON, PDB, and mmCIF parsing retains the complete
source document in memory.

## MOL and SDF input

The MOL reader accepts one V2000 or V3000 record. The SDF reader applies the same MOL parser to each
record and skips SDF data fields after `M  END`; data fields do not become molecule properties.

The supported subsets import:

- concrete element symbols represented by ChargeFW's periodic table;
- Cartesian coordinates as one conformer;
- source atom order and single, double, or triple covalent bonds; and
- formal charges from V2000 `M  CHG` records or V3000 atom `CHG=` attributes.

Query and wildcard atoms are rejected. Bond order 4 is imported as a single bond. V3000 coordinate bonds
(order 9) and hydrogen bonds (order 10) are omitted because the core graph cannot represent them; one
warning per omitted order is attached to the molecule record. Other trailing V3000 atom and bond
attributes are accepted but not used.

V2000 properties other than `M  CHG` are not interpreted. One warning per ignored three-character
property code is attached to the molecule record. Repeated ignored-property and omitted-bond warnings are
coalesced within a record and retain the first record-relative line number.

V3000 source atom ID tokens are retained exactly for mapping, including leading zeros. V2000 has no
explicit atom-ID field, so its zero-based atom-row positions are the source references.

The first MOL header line supplies the record ID and, unless another name is available, the molecule
name. MOL and SDF input always produces one conformer named `input`.

## MOL2 input

The MOL2 reader consumes `@<TRIPOS>MOLECULE` records and reads their `MOLECULE`, `ATOM`, and `BOND` data.
Unrelated sections are skipped before the atom section. The declared atom and bond counts, positive unique
atom IDs, bond references, coordinates, and required sections are validated.

The element is taken from the prefix before `.` in a standard atom type, while the atom-name column is
retained. Dummy, lone-pair, wildcard, halogen-group, hetero-group, and heavy-atom-group types are rejected.
Numeric bond types 1, 2, and 3 are imported directly; numeric type 4, aromatic type `ar`, and amide type
`am` become single bonds.

MOL2 atom charges are partial charges, not formal charges. They are parsed for validity but are not used;
all imported formal charges are zero. A record containing any nonzero input partial charge receives a
`partial_charges_ignored` warning. Each record produces one conformer named `input`.

MOL2 source atom ID tokens are retained exactly for mapping after their positive-integer uniqueness and
bond-reference semantics have been validated.

## ChargeFW molecule JSON input 1.0

Molecule input JSON is an explicit interchange format and is distinct from ChargeFW result JSON. A
document follows the [molecule input JSON Schema](../../schemas/molecule-input-1.0.schema.json) and has this
shape:

```json
{
  "schema_version": "1.0",
  "molecules": [
    {
      "id": "water-1",
      "name": "water",
      "atoms": [
        {"atomic_number": 8, "formal_charge": 0},
        {"atomic_number": 1, "formal_charge": 0}
      ],
      "bonds": [
        {"atoms": [0, 1], "order": 1}
      ],
      "conformers": [
        {"id": "model-1", "coordinates": [[0, 0, 0], [0.96, 0, 0]]}
      ]
    }
  ]
}
```

`schema_version`, `molecules`, and each molecule's non-empty `atoms` array are required. Every atom
requires integer `atomic_number` and `formal_charge` members. Optional `bonds` use zero-based atom indices
and bond orders 1, 2, or 3. Optional `conformers` contain one finite three-number coordinate row per atom;
their string `id` becomes the conformer name. Molecule `id` and `name` are optional strings. Unknown
members are ignored, but atom names are intentionally not part of schema 1.0.

Array order defines molecule, atom, bond, and conformer order. Readers may select either the first
conformer or all conformers. The whole JSON document is parsed before records are returned; malformed
JSON, another schema version, incorrect types, invalid indices, and inconsistent molecule dimensions are
errors.

## PDB and mmCIF input

PDB and mmCIF are parsed through Gemmi. PDB produces one molecule record. Each mmCIF data block containing
`_atom_site.id` produces one molecule record; blocks without coordinate data are skipped. The source block
name is retained as the mmCIF record ID. Atom-site IDs must be present and unique, but may be arbitrary
strings or integers; values such as `Csite`, `001`, and integers wider than 64 bits are retained exactly.
Parsing does not change the source document or exposed mapping.

The first selected model defines elements, formal charges, names, topology, and source mapping. Selected
atoms are restored to source row order even when Gemmi groups noncontiguous rows by hierarchy; later models
are reordered independently. Later models become conformers only when they have the same selected atom
count, element, formal charge, atom name, and author/label hierarchy identity. The selected alternate
location may differ between corresponding models and is retained per site. Readers may retain the first
model or all models. Unknown elements, empty selections, absent models, and incompatible model sequences
are errors.

Structural mappings retain original atom/site IDs, model identity, selected alternate locations, entity,
insertion code, segment where available, and distinct mmCIF author and label atom/residue/chain/sequence
identities. PDB populates its known author hierarchy and leaves unavailable label identity absent.
Each structural import also retains a component-instance partition from selected residue membership. The
record-local `component_instance_id` is the entry's zero-based position in this table, not the component
name and not a separate stored field. The component name prefers a nonempty mmCIF label residue name, then
the author residue name, then Gemmi's residue name. Each atom-index list is sorted in calculation order;
table order follows the first calculation atom. These are selected structural residues, including polymer
residues, rather than graph-connected components. The partition is stable across retained conformers, and
source hierarchy strings remain unchanged alongside it. Nonstructural imports do not infer component
instances from names or connectivity.

Structural residue template matching uses that same label-first component name and label-first atom names,
falling back to author names and then the parsed Gemmi names when a source label is unavailable. This applies
to built-in templates and mmCIF `_chem_comp_bond` rows; source-address connections continue to resolve against
the parsed structure. Names are normalized only in the private selected-residue view, so molecule atom names
and retained author/label mappings remain as imported. Duplicate canonical atom names in one selected
residue are rejected after alternate-location selection.

### Structural selection

Selection is applied before conformer validation and bonding:

- `all` retains polymer, ligand, and water residues;
- `polymers-and-ligands` excludes water residues recognized by Gemmi: `HOH`, `WAT`, `H2O`, and `DOD`; and
- `polymers` retains only residues Gemmi classifies as linear polymers. Branched, non-polymer, and water
  entities are excluded; Gemmi infers unknown entity types on the private imported structure, and residues
  that remain unclassified are excluded.

Only one alternate location is imported for each atom name in a residue: the first location in source order.
The selected atom retains its original alternate-location identity in source mapping and generated output.

### Connectivity

Structural readers provide four explicit bond strategies:

- `none` imports no bonds;
- `explicit` imports PDB `CONECT` and covalent/disulfide connections as single bonds, or mmCIF
  component bonds with their orders and covalent/disulfide connections;
- `templates` applies ChargeFW's built-in residue templates and peptide/nucleotide polymer links; and
- `hybrid` combines template and explicit bonds. Where both define a bond, mmCIF component bond orders
  are used; PDB `CONECT` records only add bonds that templates do not define.

Explicit structural aromatic bonds become single bonds. Delocalized (including `1.5`), metal,
unspecified/null, and malformed component bond orders are skipped because the core graph cannot represent
them. Templates cover common amino acids, nucleotides, water, and basic sequential peptide and nucleotide
links; they are not a complete Chemical Component Dictionary. No distance-based bond perception is performed.

For standard amino-acid residue names, `templates` and `hybrid` supplement the base templates with
N-terminal `N-H1` and `N-H3` bonds; `PRO` instead adds `N-H2` and `N-H3`. Only supplied atoms are
connected. Existing `H/H2`-only patterns (`H` for `PRO`) retain their base-template behavior.
The extension preserves supplied atoms, names, order, and formal charges. Callers supply a consistently
prepared structure; these bonds are assigned by residue and atom names, independently of chain position.

## Charge output

ChargeFW result JSON retains native floating-point precision, as do native and Python result objects.
MOL2 and mmCIF files write each coordinate as the shortest decimal number that reads back as the same
value, so coordinates read from a file keep their digits, and charges rounded to four decimal places.
The charges in such a file can therefore differ from the molecule's total charge by up to 0.00005 per
atom; use result JSON when you need exact values.
Writers validate assignment order, cardinality, targets, scope, and source references where the output
representation provides them rather than silently reordering assignments. SDF remains input-only.

### MOL2

Generated MOL2 is intended primarily for small molecules. It writes one `SMALL`/`USER_CHARGES` record for
each retained conformer, preserving native atom and bond order and joining the corresponding conformer
assignment during construction. Molecule-scoped charges are repeated for each conformer; conformer-scoped
charges are written only to their selected conformer. Record names use the caller or imported record ID,
the molecule name, or a generated fallback, with a conformer suffix when needed.

Atom and bond IDs are generated one-based values. Known atom names are retained when representable;
otherwise names are generated from the element and atom position. Atom types are element symbols, bond
types are numeric single/double/triple orders, and the generated substructure fields are `1 UNL` because
ChargeFW does not infer Tripos atom types or chemical substructures. The writer declares the generated
`UNL` substructure, but its element-only atom types are deliberately generic and may not suit consumers
that require fully classified SYBYL types.
Output rejects unsuccessful results, empty molecules, and missing or non-finite coordinates before writing.
It does not preserve source MOL2 fields, typing, substructures, comments, IDs, or lexical format.

### mmCIF

mmCIF charge output follows the SB NCBR partial atomic charges dictionary version 1.1. It adds
`_sb_ncbr_partial_atomic_charges_meta` metadata and `_sb_ncbr_partial_atomic_charges` values, plus the
dictionary declaration in `_audit_conform`. Charge rows refer to `_atom_site.id`, preserving the mapping
of structural selections and conformers.

The result-owned mmCIF writer creates a fresh minimal block for every record occurrence. It
emits `_entry`, `_audit_conform`, `_atom_site`, and the charge dictionary categories without copying source
categories or inventing component topology. Known author and label atom/residue/chain identities are kept
separate for every conformer site; unavailable names fall back deterministically to generated atom names,
`UNL`, chain `A`, and entity `1`. Output site IDs and model numbers are newly generated per block, and each
charge row refers to the generated `_atom_site.id` for its site. Molecule-scoped charges cover
every conformer and conformer-scoped charges never broadcast. Output rejects missing or non-finite
coordinates, empty records, unsuccessful results, and charges outside the dictionary's inclusive
`[-5, 5]` range before serializing the document.

Python additionally provides a strict in-memory annotation API for the unchanged Gemmi document used for
import. Its mutation and validation contract is documented with
[`attach_charges()`](../interfaces/python.md#toolkit-integrations).

### ChargeFW result JSON 1.0

Result JSON is the complete machine-readable calculation record. Its top level contains:

The complete field definitions and descriptions are available in the
[calculation result JSON Schema](../../schemas/result-1.0.schema.json).

- `schema_version`, generator identity, overall status, and document diagnostics;
- input-ordered `results`, each with source identity, record status, diagnostics, optional import mapping
  and policy, and successful charge assignments; and
- optional `calculation_provenance` containing requested and effective calculation policy and execution
  metrics.

Each assignment declares its `molecule` or `conformer` scope, charge unit `e`, and a `target` containing
the zero-based molecule index and, for conformer-specific calculations, the zero-based conformer index.
Its full-precision `charges` array follows calculation atom order and includes a full-precision
`total_charge`. Failed or cancelled records omit assignments. Diagnostics have stable severity, code, and
message fields and may include zero-based molecule, atom, bond, or conformer indices and a one-based
source line number.

For imported molecules, `input.import` records the source format, connectivity summary, and record-local
import policy. Charges follow calculation atom order. For structural inputs, generated mmCIF associates
each charge with its generated `_atom_site.id`; Python's `attach_charges()` can annotate an unchanged source
mmCIF document using the retained source mapping. Failed and cancelled results keep import policy but omit
assignments. Manually constructed molecules omit `input.import`.

Caller-supplied record and atom IDs are strings or signed 64-bit integers. Result JSON preserves those
types as JSON strings or numbers. Explicit atom IDs from a manually constructed molecule appear once as
`input.atom_ids` in calculation order.

Requested provenance records method and parameter selection, classification mode, method options,
execution request, resource thresholds, and thread limit. Import policy is record-local under
`input.import.policy`, so mixed import histories require no invocation-wide fallback.
Effective provenance records the resolved method, parameter set, complete options, execution policy, and
warnings. When supplied by the application, metrics include UTC start/end timestamps through result
finalization, parsing, applicability and computation runtimes, and peak resident memory. Durations and
memory are rounded to three decimal places.

Effective provenance may also include `fixed_ions`, an array grouped by component ID and prescribed
charge. A component ID may occur in multiple entries when its instances have different charges. Each
entry has `component_id`, `charge`, and `instances`; each monatomic instance identifies
its zero-based `molecule_index` in the original input collection and its zero-based `atom_index` in
that molecule. Component IDs use the imported label residue ID, falling back to the author residue ID.
When an explicit native fixed charge has no unambiguous imported component identity, its entry uses
`component_id: null`; it is still included and grouped with other unidentified sources of the same charge.
The optional `fixed_ions` array is an additive extension to result schema 1.0. Consumers
validating result documents with fixed ions should use the updated 1.0 schema; an older strict
schema rejects the additional effective field.
