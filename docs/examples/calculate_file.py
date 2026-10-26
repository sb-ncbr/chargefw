"""Read a molecular file, calculate charges, and print or write them per atom.

Run with: python calculate_file.py molecules.sdf --format sdf --method eem
"""

import argparse
import csv
from collections.abc import Iterable, Iterator
from pathlib import Path

import chargefw

FIELDS = ("molecule", "conformer", "atom_index", "atom_id", "atom_name", "element", "charge")


def calculate_file(
    input_path: str | Path,
    *,
    input_format: str,
    method: str,
    parameter_set: str | None = None,
) -> tuple[chargefw.MoleculeCollection, chargefw.CalculationResult]:
    """Calculate charges for every molecule in one supported input file.

    One method, parameter set, and execution mode apply to the whole collection. The parameter set
    and execution mode are selected automatically when omitted.
    """

    molecules = chargefw.io.read(input_path, format=input_format)
    result = chargefw.calculate(
        molecules,
        method=method,
        parameter_set=parameter_set,
    )
    return molecules, result


def charge_rows(
    molecules: chargefw.MoleculeCollection, result: chargefw.CalculationResult
) -> Iterator[dict[str, object]]:
    """Yield one row per atom and assignment in source order."""

    for molecule, assignments in zip(molecules, result.assignments_by_molecule):
        elements = [chargefw.periodic_table[number].symbol for number in molecule.atomic_numbers]
        for assignment in assignments:
            for atom_index, charge in enumerate(assignment.values):
                yield {
                    "molecule": assignment.molecule_index,
                    # None for geometry-independent methods, which assign charges to the molecule.
                    "conformer": assignment.conformer_index,
                    "atom_index": atom_index,
                    "atom_id": molecule.atom_ids[atom_index],
                    "atom_name": molecule.atom_names[atom_index] if molecule.atom_names else "",
                    "element": elements[atom_index],
                    "charge": float(charge),
                }


def print_charges(rows: Iterable[dict[str, object]]) -> None:
    """Print rows as an aligned table with charges rounded for reading."""

    print("molecule  conformer  atom_index  atom_id  atom_name  element      charge")
    for row in rows:
        conformer = "-" if row["conformer"] is None else row["conformer"]
        print(
            f"{row['molecule']!s:>8}  {conformer!s:>9}  {row['atom_index']!s:>10}  "
            f"{row['atom_id']!s:>7}  {row['atom_name'] or '-'!s:<9}  {row['element']!s:<7}  "
            f"{row['charge']:>10.6f}"
        )


def write_csv(path: str | Path, rows: Iterable[dict[str, object]]) -> None:
    """Write rows as CSV with full-precision charges; a missing conformer is an empty field."""

    with Path(path).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--format", required=True, choices=chargefw.io.INPUT_FORMATS)
    parser.add_argument("--method", required=True)
    parser.add_argument("--parameter-set")
    parser.add_argument(
        "--csv", type=Path, help="Write per-atom charges as CSV instead of printing"
    )
    parser.add_argument("--result-json", type=Path, help="Also write complete result JSON")
    arguments = parser.parse_args()

    molecules, result = calculate_file(
        arguments.input,
        input_format=arguments.format,
        method=arguments.method,
        parameter_set=arguments.parameter_set,
    )
    rows = charge_rows(molecules, result)
    if arguments.csv is not None:
        write_csv(arguments.csv, rows)
    else:
        print_charges(rows)
    if arguments.result_json is not None:
        chargefw.io.write(arguments.result_json, result, format="result-json")


if __name__ == "__main__":
    main()
