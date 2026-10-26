"""Calculate charges for RDKit molecules and write them as SDF atom properties.

Run with: python charge_rdkit_sdf.py input.sdf charged.sdf --method eem
"""

import argparse
import sys
from pathlib import Path

import chargefw
from chargefw.io import rdkit as chargefw_rdkit
from rdkit import Chem

PROPERTY_NAME = "ChargeFWPartialCharge"


def read_molecules(input_path: str | Path) -> list[Chem.Mol]:
    """Read SDF records with RDKit, skipping and reporting records RDKit cannot read."""

    molecules = []
    # Keep the file's explicit hydrogens; ChargeFW does not add or remove them.
    for index, molecule in enumerate(Chem.SDMolSupplier(str(input_path), removeHs=False)):
        if molecule is None:
            print(f"record {index}: skipped, RDKit could not read it", file=sys.stderr)
            continue
        molecules.append(molecule)
    if not molecules:
        raise ValueError("no SDF record could be read")
    return molecules


def convert_molecules(molecules: list[Chem.Mol]) -> chargefw.MoleculeCollection:
    """Copy RDKit molecules into one ChargeFW collection in the same order."""

    converted = []
    for index, molecule in enumerate(molecules):
        try:
            # ChargeFW needs explicit single/double/triple bonds; kekulize aromatic rings instead
            # of converting them with bond_conversion="single", which would lose the double bonds.
            Chem.Kekulize(molecule, clearAromaticFlags=True)
            name = molecule.GetProp("_Name") if molecule.HasProp("_Name") else ""
            converted.append(chargefw_rdkit.from_mol(molecule, source_name=name))
        except ValueError as error:
            raise ValueError(f"molecule {index}: {error}") from error
    return chargefw.MoleculeCollection(converted)


def attach_charges(molecules: list[Chem.Mol], result: chargefw.CalculationResult) -> None:
    """Attach each molecule's charges and the calculation model as RDKit properties."""

    plan = result.plan
    assert plan is not None, "a successful calculation reports its executed plan"
    for molecule_index, molecule in enumerate(molecules):
        chargefw_rdkit.attach_charges(
            molecule, result, molecule_index=molecule_index, property_name=PROPERTY_NAME
        )
        molecule.SetProp("ChargeFWMethod", plan.method.id)
        if plan.parameter_set is not None:
            molecule.SetProp("ChargeFWParameterSet", plan.parameter_set.id)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--method", default="eem")
    parser.add_argument("--parameter-set")
    arguments = parser.parse_args()

    molecules = read_molecules(arguments.input)
    # One collection calculation selects a single method, parameter set, and execution mode for
    # every molecule; a molecule the model cannot handle makes the calculation fail.
    result = chargefw.calculate(
        convert_molecules(molecules),
        method=arguments.method,
        parameter_set=arguments.parameter_set,
    )
    attach_charges(molecules, result)
    with Chem.SDWriter(str(arguments.output)) as writer:
        for molecule in molecules:
            writer.write(molecule)

    plan = result.plan
    assert plan is not None, "a successful calculation reports its executed plan"
    parameter_set = plan.parameter_set.id if plan.parameter_set is not None else "-"
    print(f"Model: {plan.method.id}/{parameter_set}")
    for molecule in molecules:
        charges = [atom.GetDoubleProp(PROPERTY_NAME) for atom in molecule.GetAtoms()]
        most_positive = molecule.GetAtomWithIdx(max(range(len(charges)), key=charges.__getitem__))
        print(
            f"{molecule.GetProp('_Name')!r}: total={sum(charges):+.3f} "
            f"most positive={most_positive.GetSymbol()}{most_positive.GetIdx() + 1} "
            f"({most_positive.GetDoubleProp(PROPERTY_NAME):+.3f})"
        )


if __name__ == "__main__":
    main()
