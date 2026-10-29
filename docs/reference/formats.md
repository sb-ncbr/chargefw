# File formats

This page describes what ChargeFW reads from each input format and what it writes, the same for the
[command line](../interfaces/cli.md), [Python](../interfaces/python.md), and [C++](../interfaces/cpp.md).

A file is read into molecules made of atoms in file order, bonds, formal charges, and zero or more
conformers with Cartesian coordinates. The readers change only what is described below for each
format; [Preparing input](../guide/input.md) describes what a molecule needs before a calculation. The
order of molecules and atoms is kept in the charges and in every output.

## Format overview

| Format | Accepted input | Available output | Input record model |
| --- | --- | --- | --- |
| MOL V2000/V3000 | Yes | No | One molecule, one conformer |
| SDF V2000/V3000 | Yes | No | One molecule per record, one conformer each |
| Tripos MOL2 | Yes | Yes, written new | One molecule per `MOLECULE` record, one conformer each |
| ChargeFW molecule JSON 1.0 | Yes | No | `molecules` array, zero or more conformers each |
| PDB | Yes, through Gemmi | No | One molecule, models become conformers |
| mmCIF | Yes, through Gemmi | Yes, written new | One molecule per coordinate-bearing block |
| PQR | No | Yes, for PDB and mmCIF input | Not applicable |
| ChargeFW result JSON 1.0 | No | Yes | Not applicable |

Input and output capabilities are independent. A result calculated from any supported molecular input can
be written as result JSON or, when its data is representable, as generated MOL2 or generated mmCIF. A
structure read from PDB or mmCIF can also be written as PQR. ChargeFW does not write MOL, SDF, molecule
JSON, or PDB.

Readers report malformed or unsupported data as errors. MOL, SDF, and MOL2 files may use LF or CRLF line
endings. Every record read keeps its source name, zero-based record index, record ID from the file, and
warnings, together with a mapping to the file: each atom's position and ID where the format has IDs, the
file sites of each conformer, the reading options, and whether the file had no bond data, an empty bond
list, or bonds. SDF and MOL2 files are read one record at a time; molecule JSON, PDB, and mmCIF files are
read whole into memory.

## MOL and SDF input

The MOL reader accepts one V2000 or V3000 record. The SDF reader applies the same MOL parser to each
record and skips SDF data fields after `M  END`; data fields do not become molecule properties.

The supported subsets import:

- concrete element symbols represented by ChargeFW's periodic table;
- Cartesian coordinates as one conformer;
- source atom order and single, double, or triple covalent bonds; and
- formal charges from V2000 `M  CHG` records or V3000 atom `CHG=` attributes.

Query and wildcard atoms are rejected. Bond order 4 is imported as a single bond. V3000 coordinate bonds
(order 9) and hydrogen bonds (order 10) are left out, because ChargeFW's molecules have only single,
double, and triple bonds; one warning per left-out order is attached to the record. Other trailing V3000 atom and bond
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

MOL2 atom IDs are kept exactly as written, after checking that they are unique positive integers and that
the bonds refer to them.

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
Parsing does not change the source document or exposed mapping. Partial charges in a block, such as those
of an earlier ChargeFW or ACC calculation in `_sb_ncbr_partial_atomic_charges`, are not used; the record
receives a `partial_charges_ignored` warning.

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
the author residue name, then Gemmi's residue name. Each atom-index list follows the atom order of the molecule;
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
  component bonds and covalent/disulfide `_struct_conn` connections with their orders;
- `templates` applies ChargeFW's built-in residue templates and peptide/nucleotide polymer links; and
- `hybrid` combines template and explicit bonds. Where both define a bond, mmCIF bond orders are used;
  PDB `CONECT` records only add bonds that templates do not define.

An mmCIF connection takes its order from `_struct_conn.pdbx_value_order`: `doub` and `trip` give
double and triple bonds, and any other value, or a missing one, gives a single bond.

Crystal structure files can describe bonds to symmetry-related copies of the listed atoms: a connection
whose partners have different symmetry codes (`_struct_conn.ptnr1_symmetry` and `ptnr2_symmetry` in
mmCIF, the symmetry columns of `LINK` and `SSBOND` in PDB). ChargeFW does not build those copies, so
`explicit` and `hybrid` leave such connections out; the calculation covers only the atoms in the file.
When both partner atoms are among the selected atoms, which would otherwise have been bonded wrongly,
the record receives a `symmetry_connections_ignored` warning.

Explicit structural aromatic bonds become single bonds. Delocalized (including `1.5`), metal,
unspecified/null, and malformed component bond orders are skipped, because ChargeFW's molecules have only
single, double, and triple bonds. Templates cover common amino acids, nucleotides, water, and basic sequential peptide and nucleotide
links; they are not a complete Chemical Component Dictionary. No distance-based bond perception is performed.

For standard amino-acid residue names, `templates` and `hybrid` supplement the base templates with
N-terminal `N-H1` and `N-H3` bonds; `PRO` instead adds `N-H2` and `N-H3`. Only supplied atoms are
connected. Existing `H/H2`-only patterns (`H` for `PRO`) retain their base-template behavior.
These bonds do not change atoms, names, order, or formal charges. They are assigned by residue and atom
names, independently of the position in the chain, so the structure must be prepared consistently.

## Charge output

Result JSON keeps charges at full precision, as do the Python and C++ results.
MOL2 and mmCIF files write each coordinate as the shortest decimal number that reads back as the same
value, so coordinates read from a file keep their digits, and charges rounded to four decimal places.
The charges in such a file can therefore differ from the molecule's total charge by up to 0.00005 per
atom; use result JSON when you need exact values.
Writers check that the assignments match the molecules, in order, number, targets, and scope, and refuse
to write rather than reorder them. SDF is only read.

### MOL2

Generated MOL2 is intended primarily for small molecules. It writes one `SMALL`/`USER_CHARGES` record for
each conformer, in the atom and bond order of the molecule, with that conformer's charges. Molecule-scoped charges are repeated for each conformer; conformer-scoped
charges are written only to their selected conformer. Record names use the record ID given or read, the
molecule name, or a generated name, with a conformer suffix when needed.

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

Each charge set has one metadata row, with ID `1` for the first set, and its charges cover every
conformer under that ID. In Python, one file can hold several charge sets: pass a list of results to
`chargefw.io.dumps()`, `chargefw.io.write()`, or the Gemmi functions, and the sets get IDs `1`, `2`, and
so on. Calculate all sets from the same molecules, for example one collection read once. ChargeFW
compares the input records exactly, including names, source, and import data; molecules read or built
separately are accepted only when all of these match, not merely when they are chemically identical.
Two sets with the same method and parameter set are refused, because the file does not record method
options or the execution mode that could tell them apart.

The mmCIF writer creates a new, minimal data block for every record. It
emits `_entry`, `_audit_conform`, `_atom_site`, the bond categories described below, and the charge
dictionary categories without copying source categories. Known author and label atom/residue/chain
identities are kept separate for every conformer site; unavailable names fall back deterministically to
generated atom names, `UNL`, chain `A`, and entity `1`. Output site IDs and model numbers are newly
generated per block, and each charge row refers to the generated `_atom_site.id` for its site.
Molecule-scoped charges cover every conformer and conformer-scoped charges never broadcast. Output
rejects missing or non-finite coordinates, empty records, unsuccessful results, and charges outside the
dictionary's inclusive `[-5, 5]` range before serializing the document.

The bonds used in the calculation are written with their orders (`sing`, `doub`, or `trip`), which
viewers that read mmCIF bond data can display. Bonds within a residue go to `_chem_comp_bond`, named by residue
and atom name, with the residue names listed in `_chem_comp`. Bonds between residues, such as peptide
links, go to `_struct_conn`: a bond between two sulfur atoms as `disulf`, any other bond as `covale`.
A record from SDF, MOL2, molecule JSON, or Python arrays is one residue, `UNL`, so all its bonds are in
`_chem_comp_bond`. A molecule without bonds gets no bond categories.

A `_chem_comp_bond` row applies to every residue with that name, so ChargeFW refuses to write a record
when its bonds cannot be stated that way:

- an atom name occurs twice in one residue;
- two residues with the same name have the same bond with different orders; or
- one residue has a bond that another residue with the same name and the same two atom names lacks.

A bond between residues is also refused when a partner residue has neither an author nor a label
residue number, or when a partner atom has different author and label atom or residue names; mmCIF
readers such as Gemmi would not find that atom. Nothing is written when ChargeFW refuses a record.

In Python, [`attach_charges()`](../interfaces/python.md#gemmi) instead adds the charges to the original
mmCIF document that was read, after checking that it still matches.

### PQR

PQR is the input of Poisson–Boltzmann solvers such as APBS. Each atom line has the PDB atom name,
residue name, chain, residue number, and insertion code, the coordinates in ångströms with three
decimals, the charge with four decimals, and a radius. Fields are separated by at least one space, so
programs that split lines on whitespace, as APBS does, read every line correctly; the columns therefore
differ slightly from PDB:

```text
ATOM      1  N   GLY A   12A   11.104    6.134   -6.504 -0.2500 1.5500
```

The labels are the author labels of the file, or the `label_*` values of mmCIF where author labels are
missing. Atom serial numbers count from 1. Two `REMARK` lines record the method and
parameter set and the radii.

The radius is the van der Waals radius of the atom's element, for example 1.20 Å for H, 1.70 Å for C,
1.55 Å for N, 1.52 Å for O, and 1.80 Å for S. These are not the force-field radii that PDB2PQR assigns
by atom type. ChargeFW refuses to write PQR when an element has no van der Waals radius in its periodic
table, such as Mo, W, Ru, or the lanthanides.

PQR holds one structure with one set of coordinates. ChargeFW writes it only for a result with one
molecule read from PDB or mmCIF, with one conformer. APBS ignores `MODEL` records and would treat all
models as one system, so for a file with several models, read only the first: `--conformers first` on the
command line, `conformers="first"` in Python.

APBS reads each line as words separated by whitespace. ChargeFW refuses to write labels that it would
misread:

- a chain ID that is not a single letter or other non-numeric character, such as `1`, `AA`, or `A1`;
- an atom name, residue name, or chain ID that contains whitespace, a control character, `#`, or `%`;
- an insertion code that is not a single letter.

ChargeFW checks every atom before writing, so a refused structure leaves no PQR output.

### ChargeFW result JSON 1.0

Result JSON is the complete machine-readable record of a calculation. The
[result JSON Schema](../../schemas/result-1.0.schema.json) defines every field. Its top level contains:

- `schema_version`, generator identity, overall status, and document diagnostics;
- input-ordered `results`, each with source identity, record status, diagnostics, optional import mapping
  and policy, and successful charge assignments; and
- optional `calculation_provenance` containing requested and effective calculation policy and execution
  metrics.

Each assignment declares its `molecule` or `conformer` scope, charge unit `e`, and a `target` containing
the zero-based molecule index and, for conformer-specific calculations, the zero-based conformer index.
Its `charges` array follows the atom order of the molecule, at full precision, together with the
`total_charge`. Failed or cancelled records omit assignments. Diagnostics have stable severity, code, and
message fields and may include zero-based molecule, atom, bond, or conformer indices and a one-based
source line number.

For imported molecules, `input.import` records the source format, connectivity summary, and record-local
import policy. Charges follow the atom order of the molecule. For structural inputs, generated mmCIF associates
each charge with its generated `_atom_site.id`; Python's `attach_charges()` can annotate an unchanged source
mmCIF document using the retained source mapping. Failed and cancelled results keep import policy but omit
assignments. Manually constructed molecules omit `input.import`.

Record and atom IDs you supply are strings or signed 64-bit integers, and result JSON keeps them as JSON
strings or numbers. Atom IDs of a molecule built in code appear once, as `input.atom_ids`, in atom
order.

Requested provenance records method and parameter selection, classification mode, method options,
execution request, size limits, and thread limit. The reading options are recorded per record, under
`input.import.policy`.
Effective provenance records the resolved method, parameter set, complete options, execution policy, and
warnings. When supplied by the application, metrics include UTC start/end timestamps through result
finalization, parsing, applicability and computation runtimes, and peak resident memory. Durations and
memory are rounded to three decimal places.

Effective provenance may also include `fixed_ions`, an array grouped by component ID and prescribed
charge. A component ID may occur in multiple entries when its instances have different charges. Each
entry has `component_id`, `charge`, and `instances`; each monatomic instance identifies
its zero-based `molecule_index` in the original input collection and its zero-based `atom_index` in
that molecule. Component IDs use the imported label residue ID, falling back to the author residue ID.
When a fixed charge given in C++ has no unambiguous residue name, its entry uses
`component_id: null`; it is still included and grouped with other unidentified sources of the same charge.
The optional `fixed_ions` array is an additive extension to result schema 1.0. Consumers
validating result documents with fixed ions should use the updated 1.0 schema; an older strict
schema rejects the additional effective field.
