"""Calculate point-charge dipole moments of every molecule and conformer.

Run with: python calculate_dipole_moments.py molecules.sdf --format sdf --method eem
"""

import argparse
from pathlib import Path

import chargefw
import numpy as np
from numpy.typing import NDArray

# One elementary charge times one angstrom in debye.
DEBYE_PER_E_ANGSTROM = 4.803204


def dipole_moment(
    charges: NDArray[np.float64],
    coordinates: NDArray[np.float64],
    atomic_numbers: NDArray[np.int64],
) -> NDArray[np.float64]:
    """Return the point-charge dipole vector in debye.

    Positions are taken relative to the center of nuclear charge. The dipole of a neutral molecule
    does not depend on this origin; for a charged molecule the origin is part of the result.
    """

    origin = atomic_numbers @ coordinates / atomic_numbers.sum()
    dipole: NDArray[np.float64] = DEBYE_PER_E_ANGSTROM * (charges @ (coordinates - origin))
    return dipole


def calculate_dipole_moments(
    input_path: str | Path,
    *,
    input_format: str,
    method: str,
    parameter_set: str | None = None,
) -> chargefw.CalculationResult:
    """Read every conformer, calculate charges with one model, and print dipole moments."""

    molecules = chargefw.io.read(input_path, format=input_format, conformers="all")
    if any(not molecule.has_coordinates for molecule in molecules):
        raise ValueError("every molecule needs coordinates to define a dipole moment")
    result = chargefw.calculate(molecules, method=method, parameter_set=parameter_set)

    plan = result.plan
    assert plan is not None, "a successful calculation reports its executed plan"
    parameter_set_id = plan.parameter_set.id if plan.parameter_set is not None else "-"
    print(f"Model: {plan.method.id}/{parameter_set_id}")
    print("molecule  conformer  net charge   |mu| (D)     mu_x      mu_y      mu_z")
    for molecule, assignments in zip(molecules, result.assignments_by_molecule):
        for assignment in assignments:
            # A geometry-independent method assigns one charge vector to every conformer.
            conformers = (
                range(molecule.conformer_count)
                if assignment.conformer_index is None
                else [assignment.conformer_index]
            )
            for conformer in conformers:
                dipole = dipole_moment(
                    assignment.values,
                    molecule.coordinates[conformer],
                    molecule.atomic_numbers,
                )
                print(
                    f"{assignment.molecule_index:>8}  {conformer:>9}  "
                    f"{assignment.values.sum():>+10.3f}  {np.linalg.norm(dipole):>9.3f}  "
                    f"{dipole[0]:>+8.3f}  {dipole[1]:>+8.3f}  {dipole[2]:>+8.3f}"
                )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--format", required=True, choices=chargefw.io.INPUT_FORMATS)
    parser.add_argument("--method", required=True)
    parser.add_argument("--parameter-set")
    arguments = parser.parse_args()
    calculate_dipole_moments(
        arguments.input,
        input_format=arguments.format,
        method=arguments.method,
        parameter_set=arguments.parameter_set,
    )


if __name__ == "__main__":
    main()
