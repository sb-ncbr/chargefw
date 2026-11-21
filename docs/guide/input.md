# Preparing input

ChargeFW calculates charges for the molecule exactly as you give it. Before the calculation, prepare
each molecule in the state you want charges for:

| What | Why | Where ChargeFW reads it |
| --- | --- | --- |
| All atoms, including every hydrogen | Parameter sets were fitted to molecules with explicit hydrogens | Atom records of the file |
| Bonds with their orders | Methods follow the bonds, and many parameter sets distinguish atoms by their highest bond order | SDF and MOL2 bond blocks; for PDB and mmCIF, residue templates and the bonds in the file |
| Formal charges | They set the total charge of the molecule and where charge is located | SDF `M  CHG` lines, MOL2 `UNITY_ATOM_ATTR` sections, the PDB charge columns, mmCIF `pdbx_formal_charge`, molecule JSON |
| 3D coordinates | Geometry-dependent methods need them | Atom coordinates; every model or conformer gets its own charges |

A MOL2 file states formal charges only in its optional `UNITY_ATOM_ATTR` section, which Open Babel
writes for every charged atom. ChargeFW reads an atom without such an entry as neutral; the charge
column of the atom records holds partial charges, which ChargeFW does not use. The
[format reference](../reference/formats.md) describes exactly what each reader imports.

A molecule that contains carbon but no hydrogen atoms gets a `hydrogens_missing` warning, because its
hydrogens were most likely not added: a crystal structure before protonation, or an RDKit molecule
without `Chem.AddHs()`. ChargeFW still calculates its charges. Molecules such as CO₂ or CCl₄ get the
warning as well; for them, ignore it. A molecule that lacks only some of its hydrogens gets no warning.

Many SDF files, for example the default download from PubChem, contain a 2D layout: coordinates that
only place the atoms for drawing the structural formula on a page. Their bond lengths and angles are
not those of the molecule, so ChargeFW reads an SDF or MOL record marked 2D without coordinates and
warns with `coordinates_2d_ignored`. Methods that need no geometry still apply to it;
for the others, download or generate 3D coordinates. A flat 3D molecule, such as water with all z = 0,
is not affected.

Methods that use coordinates calculate charges from atom distances, so atoms that overlap give
charges that can be off by tens of elementary charges. This happens in deposited structures with
errors and when RDKit embeds the ions of a salt on top of each other. When the method uses
coordinates, a conformer gets an `atoms_too_close` warning if two of its atoms are closer than 0.7
times the sum of their covalent radii: for example, two carbon atoms closer than 1.06 Å, or carbon
and hydrogen closer than 0.75 Å. Real bonds, including triple bonds, are longer. The warning says how
many pairs there are and names the pair that is closest relative to the sum of its radii, which is
not always the pair with the shortest distance; check the coordinates of those atoms. Its
`atom_indices`, in Python and in result JSON, lists every atom of a close pair, so that a viewer can
select them. These methods
do not calculate a conformer in which two atoms are at identical positions, unless both are
[fixed ions](fixed-ions.md); two fixed ions at one position get the warning instead.

## Small molecules

Use RDKit, or another cheminformatics toolkit, to build the molecule. Write the SMILES or structure in the
protonation state you need, add the hydrogens, and generate coordinates if the method uses geometry:

```python
from rdkit import Chem
from rdkit.Chem import AllChem

import chargefw
import chargefw.io.rdkit

mol = Chem.AddHs(Chem.MolFromSmiles("CC(=O)[O-]"))  # acetate
AllChem.EmbedMolecule(mol, randomSeed=42)
Chem.Kekulize(mol, clearAromaticFlags=True)

result = chargefw.calculate(chargefw.io.rdkit.from_mol(mol), method="eem")
print(result.assignments[0].values)
```

The charges sum to −1, the formal charge of acetate. `from_mol()` accepts single, double, and triple
bonds, so `Chem.Kekulize()` first replaces any aromatic bonds with them. The RDKit integration needs
`pip install --pre "chargefw[rdkit]"`; see [Python](../interfaces/python.md#rdkit) for
the conversion options and for writing charges back to RDKit molecules.

Files that already contain hydrogens and formal charges, such as SDF files from a docking program or a
quantum chemistry workflow, can be read directly.

## Proteins and nucleic acids

Structures from the PDB usually lack hydrogens, and the protonation of a residue depends on the pH and
on its environment. Add hydrogens and assign protonation states with a tool built for that before you
calculate charges. [Hydride](https://hydride.biotite-python.org/) writes mmCIF that ChargeFW reads
directly: it records the formal charges of the charged atoms and the bonds of the hydrogens it adds.
Ligands and non-standard residues may need separate care.

The prepared file must contain the formal charges of the charged atoms. Not every tool writes them, so
check the file before the calculation, as described [below](#check-the-molecule-before-calculating).
The PDB output of [PDB2PQR](https://server.poissonboltzmann.org/), for example, records no formal
charges and names residues after force-field conventions, such as `HID` for histidine, so ChargeFW
cannot use it as it is.

ChargeFW then takes bonds from its residue templates for the 20 standard amino acids, the nucleotides
A, C, G, U, DA, DC, DG, and DT, and water, under their PDB residue names, together with the bonds
recorded in the file. Other residues, such as ligands and modified residues, need their bonds from
the file or from a [component dictionary](#bonds-of-ligands-and-modified-residues). A file records
bonds as mmCIF `_chem_comp_bond` rows or PDB `CONECT` records. Files written by Hydride contain them.
A downloaded file may not: 1EHZ, [below](#bonds-of-ligands-and-modified-residues), has none for its
modified nucleotides. Check the bond count of your file as described in
[Check the molecule before calculating](#check-the-molecule-before-calculating).

When bond orders matter for a ligand or modified residue, use an mmCIF file that has them in
`_chem_comp_bond.value_order`, or a component dictionary. If PDB input gives a
[`bond_orders_unknown` warning](../reference/formats.md#connectivity), converting that file to
mmCIF does not recover the missing orders.

When an atom of a residue with two or more atoms gets no bond at all, ChargeFW warns with
`unbonded_residue_atoms` and names the first such residues; supply the missing bonds, in the file or
with a component dictionary, before you use a method or parameter set that depends on bonds.

Three options, named the same in Python and on the command line, control how a PDB or mmCIF file
is read: `bonds` sets where bonds come from, `selection` leaves out water or ligands, and
`conformers` reads all models or only the first. The
[format reference](../reference/formats.md#pdb-and-mmcif-input) describes them.

## Bonds of ligands and modified residues

The wwPDB [Chemical Component Dictionary](https://www.wwpdb.org/data/ccd) describes every residue and
ligand of the PDB under its residue name, with its atoms and the bonds between them. The
`chargefw-ccd` package holds these bonds in one file, the component dictionary. Install it with the
`ccd` extra and name the file when you read a structure.

The transfer RNA [1EHZ](https://www.rcsb.org/structure/1EHZ) has 11 kinds of modified nucleotides.
As downloaded from the RCSB PDB, the file has no bonds within them:

```console
$ pip install --pre "chargefw[ccd]"
$ chargefw inspect 1EHZ.cif
records: 1
record 0 (1EHZ) atoms=1821 bonds=1501 conformers=1 coordinates=true formal_charge=0 elements=C:746,N:294,O:696,Mg:6,P:76,Mn:3
$ chargefw inspect --components "$(python -c 'import chargefw_ccd; print(chargefw_ccd.path())')" 1EHZ.cif
records: 1
record 0 (1EHZ) atoms=1821 bonds=1845 conformers=1 coordinates=true formal_charge=0 elements=C:746,N:294,O:696,Mg:6,P:76,Mn:3
```

In Python, pass the file as `components`:

```python
import chargefw
import chargefw_ccd

molecules = chargefw.io.read("1EHZ.cif", format="mmcif", components=chargefw_ccd.path())
print(molecules[0].bond_count)
# 1845
```

ChargeFW uses the dictionary only when you name it. Installing `chargefw-ccd` does not change what
a script that passes no `components` reads. A script that passes `chargefw_ccd.path()` uses the
version of the dictionary that is installed; to reproduce a result, install the version the result
records, for example with `pip install chargefw-ccd==2026.11.13`.

The dictionary adds to the other sources of bonds and does not replace them:

- A residue with a template keeps the template's bonds, and a bond the file records with an order
  keeps that order.
- The dictionary gives the bonds within a residue that neither has. For a PDB file, it also gives the
  order of a `CONECT` bond whose order is unknown, when its entry has that bond.
- It does not link residues. In 1EHZ, `_struct_conn` rows in the file supply the bonds between the
  modified nucleotides and their neighbors.

The dictionary matches atoms by name. In a ligand or modified residue, hydrogens that a preparation
tool named differently from the dictionary get no bonds from it, and ChargeFW warns with
`unknown_component_atoms`. Each entry
describes one form of the component, usually the neutral one, and ChargeFW uses its bond orders as
they are. Formal charges still come from your file.

A dictionary's version is the date of the newest component in it, such as `2026.11.13`. When the
dictionary gave a molecule a bond or the order of one, the result records that version. The
[format reference](../reference/formats.md#connectivity) has the rules for each `bonds` setting.

## Check the molecule before calculating

`chargefw inspect` reports what ChargeFW read from a file. Here it reads the NMR structure of the
Trp-cage miniprotein, [1L2Y](https://www.rcsb.org/structure/1L2Y), as downloaded from the PDB:

```console
$ chargefw inspect 1L2Y.cif
records: 1
record 0 (1L2Y) atoms=304 bonds=310 conformers=38 coordinates=true formal_charge=0 elements=H:150,C:98,N:27,O:29
```

Compare the atom count, the bond count, and the formal charge with what you expect: the formal charge is
the net charge of the molecule at your pH. A different formal charge, or missing bonds, means the
structure was not prepared completely. ChargeFW calculates with the bonds as it read them and cannot
tell that one is missing, such as a disulfide bridge the file does not record, so the bond count is
worth checking.

Here the formal charge is wrong. The file contains every hydrogen, and the hydrogens describe a peptide
with charge +1: the N-terminus, Lys 8, and Arg 16 are protonated, and Asp 9 and the C-terminus are not.
The file records no formal charges, however, so ChargeFW reads every atom as neutral. NMR structures
from the PDB often look like this. ChargeFW does not derive formal charges from the hydrogens.

{download}`1l2y-h.cif` is the first model of 1L2Y prepared with Hydride at pH 7, which also writes the
formal charges:

```console
$ chargefw inspect 1l2y-h.cif
records: 1
record 0 (1L2Y) atoms=304 bonds=310 conformers=1 coordinates=true formal_charge=1 elements=H:150,C:98,N:27,O:29
```

How much missing formal charges change the charges depends on the method. EEM and SQE+qp use only the
total charge: for 1L2Y read as neutral, the charges differ by 0.003 e per atom on average. SQE+q0 also
uses which atoms carry the formal charges: read as neutral, Lys 8 NZ gets −1.08 instead of +0.14, and
Asp 9 OD2 −0.22 instead of −1.23.
