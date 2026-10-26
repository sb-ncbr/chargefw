"""Calculate charges for a structure while treating selected ions as fixed sources.

Run with: python calculate_with_fixed_ions.py structure.cif --format mmcif --ion MG
"""

import argparse
from collections.abc import Sequence
from pathlib import Path

import chargefw
import numpy as np


def calculate_structure(
    input_path: str | Path,
    *,
    input_format: str,
    ions: Sequence[str],
    method: str = "sqeqp",
    parameter_set: str | None = "SQEqp_Schindler2021_CCD_gen",
    selection: str = "all",
) -> chargefw.CalculationResult:
    """Read PDB/mmCIF input and calculate charges with the selected ions fixed."""

    molecules = chargefw.io.read(
        input_path,
        format=input_format,
        selection=selection,
    )
    return chargefw.calculate(
        molecules,
        method=method,
        parameter_set=parameter_set,
        fixed_ions=ions,
    )


def atom_label(molecule: chargefw.Molecule, atom_index: int) -> str:
    """Label an atom by its author chain, residue, sequence number, and atom name."""

    mapping = molecule.source_mapping
    labels = mapping.atoms[atom_index].structural_labels if mapping is not None else None
    if labels is None:
        return molecule.atom_names[atom_index] or str(molecule.atom_ids[atom_index])
    author = labels.author
    return f"{author.chain}/{author.residue} {author.sequence}/{author.atom}"


def print_ion_environments(result: chargefw.CalculationResult, radius: float) -> None:
    """Print each fixed ion, the charges of atoms within radius, and per-molecule balance."""

    plan = result.plan
    assert plan is not None and plan.fixed_ions is not None, "fixed ions must be resolved"
    sources = plan.fixed_ions.sources

    for molecule_index, (molecule, assignments) in enumerate(
        zip(result.molecules, result.assignments_by_molecule)
    ):
        fixed = [source for source in sources if source.molecule_index == molecule_index]
        if not fixed:
            continue
        for assignment in assignments:
            # Methods that support fixed ions require coordinates, so each assignment has a conformer.
            assert assignment.conformer_index is not None
            coordinates = molecule.coordinates[assignment.conformer_index]
            fixed_total = sum(source.charge for source in fixed)
            print(
                f"molecule={molecule_index} conformer={assignment.conformer_index}: "
                f"total={assignment.values.sum():+.3f} e "
                f"(active atoms {assignment.values.sum() - fixed_total:+.3f}, "
                f"fixed ions {fixed_total:+.3f})"
            )
            for source in fixed:
                print(
                    f"  fixed {atom_label(molecule, source.atom_index)} "
                    f"charge={source.charge:+.1f} e"
                )
                distances = np.linalg.norm(coordinates - coordinates[source.atom_index], axis=1)
                for atom_index in np.argsort(distances):
                    if atom_index == source.atom_index:
                        continue
                    if distances[atom_index] > radius:
                        break
                    element = chargefw.periodic_table[molecule.atomic_numbers[atom_index]].symbol
                    print(
                        f"    {atom_label(molecule, atom_index):<20} {element:<3} "
                        f"distance={distances[atom_index]:.2f} A "
                        f"charge={assignment.values[atom_index]:+.3f} e"
                    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--format", required=True, choices=("pdb", "mmcif"))
    ions = parser.add_mutually_exclusive_group(required=True)
    ions.add_argument("--ion", action="append", dest="ions", metavar="ID")
    ions.add_argument(
        "--common-ions",
        action="store_const",
        const=chargefw.COMMON_IONS,
        dest="ions",
        help=f"Fix the common ion catalog: {', '.join(chargefw.COMMON_IONS)}",
    )
    parser.add_argument("--method", default="sqeqp")
    parser.add_argument("--parameter-set", default="SQEqp_Schindler2021_CCD_gen")
    parser.add_argument(
        "--selection", choices=("all", "polymers-and-ligands", "polymers"), default="all"
    )
    parser.add_argument(
        "--radius", type=float, default=3.0, help="Neighbor radius around each ion in angstroms"
    )
    parser.add_argument("--result-json", type=Path)
    arguments = parser.parse_args()
    result = calculate_structure(
        arguments.input,
        input_format=arguments.format,
        ions=arguments.ions,
        method=arguments.method,
        parameter_set=arguments.parameter_set,
        selection=arguments.selection,
    )
    if result.plan is None or result.plan.fixed_ions is None:
        raise RuntimeError("the calculation did not resolve any selected fixed ions")

    for warning in result.warnings:
        print(f"Warning: {warning.message}")
    print_ion_environments(result, arguments.radius)
    if arguments.result_json is not None:
        chargefw.io.write(arguments.result_json, result, format="result-json")


if __name__ == "__main__":
    main()
