# File formats

This page describes what ChargeFW reads from each input format and what it writes, the same for the
[command line](../interfaces/cli.md), [Python](../interfaces/python.md), and [C++](../interfaces/cpp.md).

A file is read into molecules made of atoms in file order, bonds, formal charges, and zero or more
conformers with Cartesian coordinates. The readers change only what is described below for each
format; [Preparing input](../guide/input.md) describes what a molecule needs before a calculation. The
order of molecules and atoms is kept in the charges and in every output.

## Format overview

| Format | Read | Written | What a file holds when read |
| --- | --- | --- | --- |
| MOL V2000/V3000 | Yes | No | One molecule, one conformer (none for a 2D layout) |
| SDF V2000/V3000 | Yes | No | One molecule per record, one conformer each (none for a 2D layout) |
| Tripos MOL2 | Yes | Yes | One molecule per `MOLECULE` record, one conformer each |
| ChargeFW molecule JSON 1.0 | Yes | No | One molecule per entry of `molecules`, any number of conformers |
| PDB | Yes | No | One molecule; its models become conformers |
| mmCIF | Yes | Yes | One molecule per data block with coordinates; its models become conformers |
| PQR | No | Yes, for PDB and mmCIF input | |
| ChargeFW result JSON 1.0 | No | Yes | |

Charges calculated from any input can be written as result JSON, and as MOL2 or mmCIF when the
molecules have coordinates. Structures read from PDB or mmCIF can also be written as PQR. ChargeFW
does not write MOL, SDF, molecule JSON, or PDB.

Malformed or unsupported data in a file is an error. MOL, SDF, and MOL2 files may use LF or CRLF
line endings.

MOL, SDF, MOL2, PDB, and mmCIF files must be ASCII text throughout, including comments and SDF data
items that ChargeFW does not read. Any other character, such as an accented letter or a UTF-8 byte
order mark, is an error that names its line; remove it before reading the file. Molecule JSON files
are UTF-8.

For every molecule read, ChargeFW keeps where it came from: the file name, the
record's position in the file and its ID, the warnings, and for every atom its position and ID in
the file, so that charges can be matched back to it.

For MOL, SDF, and MOL2 input, diagnostic line numbers count lines from the start of
the input, beginning at 1, including blank lines and skipped data. They do not restart for each
record. Parsing errors also identify the record, numbered from 1, and its ID when available.
An unexpected end of input points to the next expected line.

## MOL and SDF input

A MOL file holds one V2000 or V3000 record, as defined in the BIOVIA
[CTfile Formats](https://www.3ds.com/assets/invest/2026-01/biovia_chemistry_2026_ctfile_formats.pdf)
specification; an SDF file holds any number of them. The data fields after `M  END` in SDF are
skipped. From each record, ChargeFW reads:

- the element of every atom, in file order;
- the coordinates, as one conformer;
- single, double, and triple bonds; and
- formal charges: in V3000 from the `CHG=` attribute of the atoms, in V2000 from `M  CHG` lines or,
  in a record without `M  CHG` and `M  RAD` lines, from the charge field of the atom lines.

Everything else in a record is ignored without a warning, for example isotopes, stereochemistry,
atom aliases, S-groups, and V3000 blocks other than `ATOM` and `BOND`. Because S-groups are ignored,
a polymer repeat unit is calculated as the atoms written in the file.

Molfiles have no atom names, so imported atoms have none. MOL2 and mmCIF output give them names from
the element and atom position, such as `C1`, `C2`, and `O3`.

A record marked `2D` in columns 21-22 of its second line is a 2D layout: its coordinates only place
the atoms for drawing the structural formula on a page, and its bond lengths and angles are not
those of a molecule. ChargeFW imports such a record without coordinates and attaches a
`coordinates_2d_ignored` warning. Coordinates with all z = 0 in a record not marked `2D` are kept.

Query and wildcard atoms are errors. Aromatic bonds (order 4) are read as single bonds, and the
record gets one `aromatic_bonds_converted` warning that counts them. V3000 coordinate bonds (order
9) and hydrogen bonds (order 10) are left out, because ChargeFW has only single, double, and triple
bonds; the record gets one `v3000_bond_order_ignored` warning for each of the two kinds, with the
line number of the first such bond. A V3000 record with no bonds may omit the `BOND` block or include
an empty one.

V3000 atom IDs are kept exactly as written, including leading zeros. V2000 has no atom IDs; an
atom's ID is its position in the record, counted from zero.

The first line of a record gives its ID and, unless another name is available, the molecule name. A
MOL or SDF record produces one conformer named `input`, or none for a 2D layout.

## MOL2 input

Each `@<TRIPOS>MOLECULE` record is one molecule. ChargeFW reads its `MOLECULE`, `ATOM`, `BOND`, and
`UNITY_ATOM_ATTR` sections, in any order, and skips the others. Wrong atom or bond counts, atom IDs
that are not unique positive integers, bonds to atoms that do not exist, invalid coordinates, and
missing sections are errors, as are repeated `ATOM` or `BOND` sections within one molecule. A
section marker starts with `@` in the first column; whitespace after it is ignored.

The element is the part of the atom type before the `.`, such as C in `C.ar`; the atom name is kept.
Types that do not name one element, such as dummy atoms, lone pairs, wildcards, and the groups
`Hal`, `Het`, and `Hev`, are errors. Numeric bond types 1, 2, and 3 are imported directly; numeric
type 4, aromatic type `ar`, and amide type `am` become single bonds. Aromatic bonds, types 4 and
`ar`, give the record an `aromatic_bonds_converted` warning.

Formal charges come from the `charge` attributes in `UNITY_ATOM_ATTR`, such as `charge -1` for atom 9:

```text
@<TRIPOS>UNITY_ATOM_ATTR
9 1
charge -1
```

Atoms without a `charge` attribute are neutral, and other attributes are ignored. The charge column
of the atom records holds partial charges, not formal charges; it must contain numbers but is not
used, and a record with any nonzero partial charge gets a `partial_charges_ignored` warning. Each
record produces one conformer named `input`.

MOL2 atom IDs are kept exactly as written, after checking that they are unique positive integers and that
the bonds refer to them.

## ChargeFW molecule JSON input 1.0

Molecule JSON is ChargeFW's own format for giving molecules as JSON; it is not the result JSON. A
document follows the [molecule input JSON Schema](../../schemas/molecule-input-1.0.schema.json) and
looks like this:

```json
{
  "schema_version": "1.0",
  "molecules": [
    {
      "id": "water-1",
      "name": "water",
      "atoms": [
        {"atomic_number": 8, "formal_charge": 0},
        {"atomic_number": 1, "formal_charge": 0},
        {"atomic_number": 1, "formal_charge": 0}
      ],
      "bonds": [
        {"atoms": [0, 1], "order": 1},
        {"atoms": [0, 2], "order": 1}
      ],
      "conformers": [
        {"id": "model-1", "coordinates": [[0, 0, 0], [0.96, 0, 0], [-0.24, 0.93, 0]]}
      ]
    }
  ]
}
```

`schema_version`, `molecules`, and each molecule's `atoms`, with at least one atom, are required.
Every atom needs an integer `atomic_number` and `formal_charge`. `bonds` is optional; it gives atom
indices, counted from zero, and bond orders 1, 2, or 3. `conformers` is optional; each conformer has
one row of three finite coordinates per atom, and its `id` becomes the conformer's name. A
molecule's `id` and `name` are optional strings. Atom names cannot be given in version 1.0.

Molecules, atoms, bonds, and conformers keep the order of their arrays. You can read all conformers
of each molecule or only the first. Invalid JSON, another schema version, wrong types, members not
listed here, indices out of range, and conformers whose number of rows differs from the number of
atoms are errors. A misspelled member, such as `bond` for `bonds`, therefore stops the reading
instead of giving a molecule without bonds.

## PDB and mmCIF input

ChargeFW reads PDB and mmCIF files with Gemmi. A PDB file is one molecule. In mmCIF, every data
block with atom coordinates (`_atom_site`) is one molecule, whose ID is the block name; blocks
without coordinates are skipped. Every mmCIF atom-site row needs a unique ID, which can be any
string or integer; IDs such as `Csite`, `001`, or integers with more than 64 bits are kept exactly
as written. ChargeFW does not change the file. Partial charges in it, such as those of an earlier
ChargeFW or ACC calculation in `_sb_ncbr_partial_atomic_charges`, are not used, and the record gets
a `partial_charges_ignored` warning.

The first model gives the atoms of the molecule, in the order of the rows in the file, with their
elements, formal charges, and names. Every further model must contain the same atoms: the same
elements, formal charges, atom names, residues, and chains, in the same order. Each such model
becomes another conformer, and a model that differs is an error. The alternate location chosen for
an atom may differ between models. `--conformers first` (`conformers="first"` in Python) reads only
the first model. Unknown elements and a selection that leaves no atoms are errors.

For every atom, ChargeFW keeps its ID, model, alternate location, entity, insertion code, and
segment, and its atom, residue, and chain names as the file gives them. mmCIF has two sets of names,
author names (`auth_*`) and label names (`label_*`), and both are kept; PDB has only author names.

ChargeFW also groups the atoms into residues; in Python, `molecule.source_mapping.components` lists
them, each with its name and atoms. A residue's name is its mmCIF label name where the file has one,
otherwise its author name. Residue templates and mmCIF `_chem_comp_bond` rows are matched by this
residue name and by the atom names, taken the same way. Two atoms with the same name in one residue,
after alternate locations are chosen, are an error.

### Structural selection

`selection` (`--selection` on the command line) chooses which residues are read. Bonds and further
models are then read for these atoms only.

- `all` keeps polymer, ligand, and water residues;
- `polymers-and-ligands` leaves out water: the residues `HOH`, `WAT`, `H2O`, and `DOD`; and
- `polymers` keeps only residues of polymer entities, such as protein and nucleic-acid chains, as
  the file defines its entities. Where it does not, Gemmi assigns them, and a residue in the same
  chain as a polymer, such as a ligand, can then count as part of it; residues Gemmi cannot assign
  are left out.

For each atom name in a residue, only one alternate location is read: the first one in the file. Its
alternate-location code is kept and written to the output.

### Connectivity

`bonds` (`--bonds` on the command line) chooses where bonds come from:

- `none`: no bonds;
- `explicit`: the bonds in the file. From PDB, `CONECT` records and covalent and disulfide
  connections, all as single bonds; from mmCIF, `_chem_comp_bond` rows and covalent and disulfide
  `_struct_conn` connections, with their orders;
- `templates`: ChargeFW's residue templates and the links between consecutive amino acids or
  nucleotides; and
- `hybrid`, the default: both. Where both give a bond, the mmCIF bond order is used; PDB `CONECT`
  records only add bonds that the templates do not have.

An mmCIF connection takes its order from `_struct_conn.pdbx_value_order`: `doub` and `trip` give
double and triple bonds, and any other value, or a missing one, gives a single bond.

Crystal structure files can describe bonds to symmetry-related copies of the listed atoms: a connection
whose partners have different symmetry codes (`_struct_conn.ptnr1_symmetry` and `ptnr2_symmetry` in
mmCIF, the symmetry columns of `LINK` and `SSBOND` in PDB). ChargeFW does not build those copies, so
`explicit` and `hybrid` leave such connections out; the calculation covers only the atoms in the file.
When both partner atoms are among the selected atoms, which would otherwise have been bonded wrongly,
the record receives a `symmetry_connections_ignored` warning.

With any strategy other than `none`, a record in which some residue of two or more atoms has an atom
without bonds receives one `unbonded_residue_atoms` warning. It counts these atoms and residues and
names up to three residues, such as `MSE A 4 (N, CA, C, O, ...)`. Single-atom residues, such as ions,
are not reported.

Aromatic bonds in the file become single bonds, and the record gets an `aromatic_bonds_converted`
warning. Bonds with other orders in `_chem_comp_bond`, such as delocalized (`delo`, `1.5`), metal,
or missing orders, are skipped, because ChargeFW has only single, double, and triple bonds. The
templates cover the 20 standard amino acids, the nucleotides A, C, G, U, DA, DC, DG, and DT, and
water, not the whole Chemical Component Dictionary. ChargeFW never adds bonds from the distances
between atoms.

For the standard amino acids, `templates` and `hybrid` also bond the N-terminal hydrogens `H1` and
`H3` to the nitrogen, and for proline `H1`, `H2`, and `H3`. This happens only for hydrogens present
in the file, and only in a residue that the templates do not link to the previous residue by a
peptide bond. The hydrogens `H` and `H2` (for proline `H`) are bonded by the standard templates
already.

These bonds do not change atoms, names, order, or formal charges. They are assigned by residue and atom
names and, for the links between residues and the N-terminal hydrogens, by the order of residues in
the chain; no geometry is used, so the structure must be prepared consistently.

Histidine and arginine templates follow the residue's protonation state. The templates come from the
Chemical Component Dictionary, which gives histidine a double bond `ND1=CE1` and arginine `CZ=NH2`.
A residue whose hydrogens and formal charges do not fit those bonds gets the structure from the wwPDB
protonation variants that fits: `CE1=NE2` for histidine, `CZ=NH1` or `NE=CZ` for arginine. A structure
fits when, among the histidine ring atoms or the arginine `CD`, `NE`, `CZ`, `NH1`, and `NH2`, every carbon
has bond orders summing to 4 and every nitrogen to 3 plus its formal charge. For example:

| Residue | Hydrogens | Formal charges | Bond |
| --- | --- | --- | --- |
| Histidine, neutral, protonated on ND1 | `HD1` | none | `CE1=NE2` |
| Histidine, as wwPDB variant `HIS_LL_DHE2` | `HD1` | ND1 +1, NE2 −1 | `ND1=CE1` |
| Histidine, neutral, protonated on NE2 | `HE2` | none | `ND1=CE1` |
| Histidine, cation | `HD1`, `HE2` | NE2 +1 | `CE1=NE2` |
| Arginine, neutral, without `HH12` | `HE`, `HH11`, `HH21`, `HH22` | none | `CZ=NH1` |
| Arginine, neutral, without `HE` | `HH11`, `HH12`, `HH21`, `HH22` | none | `NE=CZ` |

A missing formal charge counts as 0, so a histidine or arginine cation without its formal charge
fits no structure and keeps the template bonds; check the formal charges as described in
[Preparing input](../guide/input.md). A residue with an incomplete group, for example without `CB`
or without a hydrogen of `CD`, also keeps the template bonds. Under `hybrid`, bond orders from mmCIF
`_chem_comp_bond` rows take precedence over the template's.

## Charge output

Result JSON, MOL2, PQR, and mmCIF output, including charges attached to a Gemmi document, name
ChargeFW and the library version as the program that wrote them.

Result JSON keeps charges at full precision, as do the Python and C++ results. MOL2 and mmCIF files
write each coordinate with the fewest digits that read back as exactly the same value, and charges
rounded to four decimal places, so the charges in such a file can differ from the molecule's total
charge by up to 0.00005 per atom; use result JSON when you need exact values. ChargeFW writes a file
only when the charges belong to the molecules being written, in the same order; it never reorders
them.

### MOL2

MOL2 output is meant for small molecules. It writes one `SMALL`/`USER_CHARGES` record for each
conformer, with the atoms and bonds in the order of the molecule and the charges of that conformer.
Charges from a method that does not use coordinates are the same in every record. A record is named
after the record ID, the molecule name, or its position, with the conformer number added when there
are several.

Atoms and bonds are numbered from 1. Atom names are kept unless they contain spaces; otherwise, and
for atoms without a name, they are made from the element and the atom's position, such as `C1`. Bond
types are the numbers 1, 2, and 3, and all atoms are in one substructure, `1 UNL`.

Atom types are SYBYL types derived from the elements, bonds, bond orders, and formal charges:

| Element | Types |
| --- | --- |
| C | `C.1` (triple bond or two double bonds), `C.2` (double bond), `C.cat` (amidinium and guanidinium carbon), `C.3` |
| N | `N.1` (triple bond), `N.2` (double bond), `N.pl3` (charged or pentavalent with a double bond, as in nitro; or bonded to a carbon or nitrogen with a multiple bond, as in anilines), `N.am` (amide), `N.4` (four bonds and charge +1), `N.3` |
| O | `O.co2` (carboxylate or phosphate oxygens: one double-bonded, the others single-bonded with charge −1), `O.2` (double bond), `O.3` |
| S | `S.O2` (two oxygens bonded to sulfur only, by a double bond or as S⁺–O⁻), `S.O` (one such oxygen), `S.2` (double bond), `S.3` |
| P | `P.3` |
| Others | Element symbol, such as `H`, `Cl`, `Na`, or `Zn` |

An aromatic ring is written in the Kekulé form ChargeFW holds, with types such as `C.2` and bond orders
1 and 2; the types `C.ar` and `N.ar` and the bond type `ar` are not written. RDKit perceives the ring as
aromatic when it reads the file. Other programs may keep the Kekulé form.

Atoms with a nonzero formal charge are listed in the `UNITY_ATOM_ATTR` section with two attributes
holding the charge, `charge`, which Open Babel and ChargeFW read, and `AtomExpr`, which RDKit reads:

```text
@<TRIPOS>UNITY_ATOM_ATTR
9 2
charge -1
AtomExpr -1
```

ChargeFW does not write failed results, molecules without atoms, or molecules without valid
coordinates. The file is written new: atom types, substructures, comments, and IDs of a MOL2 input
file are not kept.

### mmCIF

mmCIF charge output follows the SB NCBR partial atomic charges dictionary version 1.1. It adds
`_sb_ncbr_partial_atomic_charges_meta` metadata and `_sb_ncbr_partial_atomic_charges` values, plus
the dictionary declaration in `_audit_conform`. Each charge refers to the `_atom_site.id` of its
atom in its model.

Each charge set has one metadata row, with ID `1` for the first set, and its charges cover every
conformer under that ID. In Python, one file can hold several charge sets: pass a list of results to
`chargefw.io.dumps()`, `chargefw.io.write()`, or the Gemmi functions, and the sets get IDs `1`, `2`, and
so on. Calculate all sets from the same molecules, for example one collection read once. ChargeFW
compares the input records exactly, including names, source, and import data; molecules read or built
separately are accepted only when all of these match, not merely when they are chemically identical.
Two sets with the same method and parameter set are refused, because the file does not record method
options or the execution mode that could tell them apart.

For every molecule, ChargeFW writes a new data block with `_entry`, `_audit_conform`, `_atom_site`,
the bond categories described below, and the charge categories; nothing else from the input file is
copied. Author and label names are written as they were read. Missing names are replaced by
generated atom names, residue `UNL`, chain `A`, and entity `1`. Atom IDs and model numbers are
numbered anew, and each charge refers to the new `_atom_site.id` of its atom. Charges from a method
that does not use coordinates are written for every model. ChargeFW does not write failed results,
molecules without atoms or without valid coordinates, or charges outside −10 to +10, the range of
the dictionary.

The bonds used in the calculation are written with their orders (`sing`, `doub`, or `trip`), which
viewers that read mmCIF bond data can display. Bonds within a residue go to `_chem_comp_bond`, named
by residue and atom name. Bonds between residues, such as peptide links, go to `_struct_conn`: a bond
between two sulfur atoms as `disulf`, any other bond as `covale`.
A record from SDF, MOL2, molecule JSON, or Python arrays is one residue, `UNL`, so all its bonds are in
`_chem_comp_bond`. A molecule without bonds gets no bond categories.

A `_chem_comp_bond` row applies to every residue with that name. Histidines or arginines in different
protonation states, such as one histidine protonated on ND1 and another on NE2, have their double bonds
on different atoms (see [Connectivity](#connectivity)). ChargeFW writes the bonds whose order differs
with the delocalized order `delo`: in histidine `ND1-CE1` and `CE1-NE2`, in arginine `NE-CZ`, `CZ-NH1`,
or `CZ-NH2`. When you read the file with the
`hybrid` bond strategy, the default, or with `templates`, ChargeFW gives each residue its bond orders
again from its hydrogens and formal charges, so the charges you calculate from the file are the same.
With `explicit`, these bonds are missing.

ChargeFW refuses to write a record when its bonds cannot be stated per residue name:

- an atom name occurs twice in one residue;
- two residues with the same name have the same bond with different orders, other than the histidine
  and arginine bonds above; or
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
- a residue number that is not an integer, such as the hybrid-36 number `A000` that PDB files use for
  residue 10000 and above (read the structure from mmCIF instead); or
- an insertion code that is not a single letter.

ChargeFW checks every atom before writing, so a refused structure leaves no PQR output.

### ChargeFW result JSON 1.0

Result JSON is the complete record of a calculation, for programs and for reporting. The
[result JSON Schema](../../schemas/result-1.0.schema.json) defines every field. At the top level it
has:

- `schema_version`, the `generator` with its name and version, the overall `status`, and the
  `diagnostics` that concern no single molecule;
- `results`, one record per input molecule in input order, each with `input`, where the molecule
  came from, its `status` and `diagnostics`, and the `assignments` of a successful calculation; and
- `calculation_provenance`: what was requested, what was used, and the run time and memory.

Later versions 1.x add members and values but do not remove or rename any, so a program that reads
1.0 can read every 1.x document. Check only that `schema_version` starts with `1.`, ignore members
you do not know, and accept values you do not know, such as a new status or input format. Two lists
stay fixed in every 1.x version: an assignment's `scope` and a record's `source_connectivity`. A
version 2.0 would be incompatible.

Each assignment has its `scope`, `molecule` or `conformer`, the `charge_unit` `e`, and a `target`
with the molecule index and, for methods that use coordinates, the conformer index, both counted
from zero. `charges` follows the atom order of the molecule, at full precision, and `total_charge`
is their sum. Failed and cancelled records have no assignments.

A diagnostic has a `severity`, a `code`, such as `atoms_too_close`, and a `message`. It can give the
molecule, atom, bond, and conformer it concerns, counted from zero, and the `line` of the file,
counted from one. A diagnostic about several atoms or bonds lists all of them in `atom_indices` or
`bond_indices`; `atom_index` and `bond_index` are the ones its message names. A viewer can select
the listed atoms and bonds to show where the problem is. For PDB and mmCIF input, the message also
names the residue of each atom it mentions, as PDB files and viewers do, for example
`atoms 3707 and 3708 are in PRO D 58`. The command line prints the same messages.

For a molecule read from a file, `input.import` records the format, the reading options (`policy`),
and whether the file contained bonds (`source_connectivity`); failed and cancelled records keep it
too. Molecules built in code have no `input.import`, and their atom IDs are listed once, in
`input.atom_ids`. Record and atom IDs are strings or signed 64-bit integers, written as JSON strings
or numbers.

`calculation_provenance.requested` records what you asked for, with the names of the arguments of
`chargefw.calculate()` in Python: `method`, `parameter_set`, `method_scope`, `method_options` (the
options given for the named method), `parameter_matching`, `execution`, `radius` in Å,
`cutoff_threshold` and `cover_threshold` in atoms or `"unlimited"`, `threads`, and the residue names
in `fixed_ions`. `effective` records the method, parameter set, all option values, and execution
that were used. Warnings about the plan, such as `resource_threshold_exceeded`, are diagnostics with
severity `warning`: in the record of the molecule they concern, or in the document's diagnostics.
`execution_metrics` gives the start and end time in UTC, the time of each step in seconds, rounded to
three decimal places, and the peak memory of the process in bytes.

With fixed ions, `effective` also has `fixed_ions`: one entry per residue name and charge, with
`component_id`, the residue name (mmCIF label name first), the `charge`, and `instances`, each with
the `molecule_index` and `atom_index` of one ion. An ion given in C++ without a residue name has
`component_id: null`.
