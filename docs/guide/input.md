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
`pip install --pre "chargefw[rdkit]"`; see [Python](../interfaces/python.md#toolkit-integrations) for
the conversion options and for writing charges back to RDKit molecules.

Files that already contain hydrogens and formal charges, such as SDF files from a docking program or a
quantum chemistry workflow, can be read directly.

## Proteins and nucleic acids

Structures from the PDB usually lack hydrogens, and the protonation of a residue depends on the pH and
on its environment. Add hydrogens and assign protonation states with a tool built for that, such as
[Hydride](https://hydride.biotite-python.org/) or [PDB2PQR](https://server.poissonboltzmann.org/),
before you calculate charges. Which tool fits depends on your structure; ligands and non-standard
residues in particular may need separate care.

The prepared file must contain the formal charges of the charged atoms. Not every tool writes them,
so check the file before the calculation, as described [below](#check-the-molecule-before-calculating).

ChargeFW then takes bonds from its residue templates for common amino acids, nucleotides, and water,
together with the bonds recorded in the file. For other residues, the file must contain their bonds:
mmCIF `_chem_comp_bond` rows or PDB `CONECT` records. The
[format reference](../reference/formats.md#connectivity) describes the bond options.

## Check the molecule before calculating

`chargefw inspect` reports what ChargeFW read from a file:

```console
$ chargefw inspect 1crn-h.cif
records: 1
record 0 (STRUCTURE) atoms=642 bonds=652 conformers=1 coordinates=true formal_charge=0 elements=H:315,C:202,N:55,O:64,S:6
```

Compare the atom count, the bond count, and the formal charge with what you expect: the formal charge is
the net charge of the molecule at your pH. A different formal charge, or missing bonds, means the
structure was not prepared completely.
