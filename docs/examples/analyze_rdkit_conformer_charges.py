"""Analyze calculated charge variation across an RDKit conformer set.

Run with: python analyze_rdkit_conformer_charges.py "CCO" --conformers 20
"""

import argparse

import chargefw
import numpy as np
from chargefw.io import rdkit as chargefw_rdkit
from rdkit import Chem
from rdkit.Chem import AllChem


def prepare_demo_ensemble(smiles: str, conformer_count: int) -> Chem.Mol:
    """Generate an unoptimized, deterministic ETKDG ensemble for demonstration."""

    if conformer_count < 2:
        raise ValueError("conformer_count must be at least 2")
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError(f"invalid SMILES: {smiles}")

    # ChargeFW's default RDKit import requires explicit single/double/triple bonds.
    # Kekulize aromatic rings to retain single/double bond orders rather than using
    # bond_conversion="single", which would turn every aromatic bond into a single bond.
    Chem.Kekulize(molecule, clearAromaticFlags=True)
    # Charge calculation needs explicit hydrogens: omitting them calculates charges for an incomplete molecule.
    # Add them before embedding so their coordinates are generated too.
    molecule = Chem.AddHs(molecule)

    # Generate 3D conformers with RDKit's ETKDGv3 distance-geometry settings:
    # https://www.rdkit.org/docs/GettingStartedInPython.html#working-with-3d-molecules
    embedding = AllChem.ETKDGv3()
    # Fix the random seed for repeatable sampling
    embedding.randomSeed = 42
    conformer_ids = AllChem.EmbedMultipleConfs(molecule, numConfs=conformer_count, params=embedding)
    if len(conformer_ids) < 2:
        raise RuntimeError(f"RDKit generated only {len(conformer_ids)} conformers")
    return molecule


def atom_charge_statistics(
    molecule: Chem.Mol, *, source_name: str = ""
) -> tuple[list[tuple[str, float, float, float, float]], chargefw.CalculationResult]:
    """Calculate per-atom mean, standard deviation, minimum, and maximum charge."""

    native_molecule = chargefw_rdkit.from_mol(molecule, source_name=source_name)
    result = chargefw.calculate(
        native_molecule,
        method="eem",
    )
    charges = np.stack([assignment.values for assignment in result.assignments])
    statistics = []
    for atom in molecule.GetAtoms():
        values = charges[:, atom.GetIdx()]
        statistics.append(
            (
                f"{atom.GetSymbol()}{atom.GetIdx() + 1}",
                float(np.mean(values)),
                float(np.std(values)),
                float(np.min(values)),
                float(np.max(values)),
            )
        )
    return statistics, result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("smiles", nargs="?", default="CCO")
    parser.add_argument("--conformers", type=int, default=20)
    parser.add_argument("--top", type=int, default=5)
    arguments = parser.parse_args()

    molecule = prepare_demo_ensemble(arguments.smiles, arguments.conformers)
    statistics, result = atom_charge_statistics(molecule, source_name=arguments.smiles)
    parameter_set = result.plan.parameter_set.id if result.plan.parameter_set else "-"
    print(f"Method: {result.plan.method.id}/{parameter_set}")
    print(f"Conformers: {len(result.assignments)} equally weighted, unoptimized ETKDGv3 structures")
    print("Atom          mean         std         min         max")
    for label, mean, standard_deviation, minimum, maximum in sorted(
        statistics, key=lambda item: item[2], reverse=True
    )[: arguments.top]:
        print(f"{label:<6} {mean:12.6f} {standard_deviation:11.6f} {minimum:11.6f} {maximum:11.6f}")


if __name__ == "__main__":
    main()
