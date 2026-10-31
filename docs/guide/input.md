# Preparing input

ChargeFW calculates charges for the molecule exactly as you give it. Before the calculation, prepare
each molecule in the state you want charges for:

| What | Why | Where ChargeFW reads it |
| --- | --- | --- |
| All atoms, including every hydrogen | Parameter sets were fitted to molecules with explicit hydrogens | Atom records of the file |
| Bonds with their orders | Methods follow the bonds, and many parameter sets distinguish atoms by their highest bond order | SDF and MOL2 bond blocks; for PDB and mmCIF, residue templates and the bonds in the file |
| Formal charges | They set the total charge of the molecule and where charge is located | SDF `M  CHG` lines, the PDB charge columns, mmCIF `pdbx_formal_charge`, molecule JSON |
| 3D coordinates | Geometry-dependent methods need them | Atom coordinates; every model or conformer gets its own charges |

MOL2 files have no field for formal charges, so ChargeFW reads every atom of a MOL2 file as neutral.
Use SDF, mmCIF, or molecule JSON for charged molecules. The
[format reference](../reference/formats.md) describes exactly what each reader imports.

Many SDF files, for example the default download from PubChem, contain a 2D layout: coordinates that
only place the atoms for drawing the structural formula on a page. Their bond lengths and angles are
not those of the molecule, so ChargeFW reads an SDF or MOL record marked 2D without coordinates and
warns with `coordinates_2d_ignored`. Methods that need no geometry still apply to it;
for the others, download or generate 3D coordinates. A flat 3D molecule, such as water with all z = 0,
is not affected.

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
recorded in the file. For other residues, the file must contain their bonds: mmCIF `_chem_comp_bond`
rows or PDB `CONECT` records. mmCIF files downloaded from the PDB and files written by Hydride contain
them. When an atom of a residue with two or more atoms gets no bond at all, ChargeFW warns with
`unbonded_residue_atoms` and names the first such residues; add the missing bonds to the file before
you use a method or parameter set that depends on bonds.

Three options control how a PDB or mmCIF file is read: `--structural-bonds` (`bonds` in Python) sets
where bonds come from, `--structural-selection` (`selection`) leaves out water or ligands, and
`--conformers` (`conformers`) reads all models or only the first. The
[format reference](../reference/formats.md#pdb-and-mmcif-input) describes them.

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
