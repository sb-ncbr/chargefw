"""Inspect the normalized molecules imported by ChargeFW.

Run with: python inspect_molecules.py input.sdf --format sdf --classifier bonded
"""

import argparse
from collections import Counter
from pathlib import Path

import chargefw

CLASSIFIERS = ("plain", "hbo", "bonded")


def atom_types(molecule: chargefw.Molecule, classifier: str) -> list[tuple[str, str]]:
    """Return the element and classifier type of every atom, following the parameter-set format."""

    symbols = [chargefw.periodic_table[number].symbol for number in molecule.atomic_numbers]
    neighbors: list[list[str]] = [[] for _ in symbols]
    highest_order = [0] * len(symbols)
    for first, second, order in molecule.bonds:
        neighbors[first].append(symbols[second])
        neighbors[second].append(symbols[first])
        highest_order[first] = max(highest_order[first], int(order))
        highest_order[second] = max(highest_order[second], int(order))

    if classifier == "plain":
        types = ["*"] * len(symbols)
    elif classifier == "hbo":
        types = [str(order) for order in highest_order]
    else:
        types = ["".join(sorted(bonded)) for bonded in neighbors]
    return list(zip(symbols, types))


def inspect_molecules(
    input_path: str | Path,
    *,
    input_format: str,
    selection: str = "all",
    bonds: str | None = None,
    conformers: str = "all",
    classifier: str | None = None,
) -> None:
    """Print source identities, dimensions, formal charges, and optional atom-type counts."""

    molecules = chargefw.io.read(
        input_path,
        format=input_format,
        selection=selection,
        bonds=bonds,
        conformers=conformers,
    )
    print(f"Collection: {molecules.name or '<unnamed>'} ({len(molecules)} molecule(s))")

    for index, molecule in enumerate(molecules):
        print(
            f"[{index}] record={molecule.record_index} id={molecule.record_id!r} "
            f"name={molecule.name!r}"
        )
        print(
            f"    atoms={molecule.atom_count} bonds={molecule.bond_count} "
            f"conformers={molecule.conformer_count} "
            f"formal_charge={int(molecule.formal_charges.sum())}"
        )
        print(f"    conformer_names={molecule.conformer_names!r}")
        if classifier is None:
            continue
        counts = Counter(atom_types(molecule, classifier))
        for (symbol, atom_type), count in sorted(counts.items()):
            # An atom without bonds has an empty bonded type; quote it to keep it visible.
            shown_type = atom_type or '""'
            print(f"    {symbol:<3} {classifier}={shown_type}: {count}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--format", choices=chargefw.io.INPUT_FORMATS, required=True)
    parser.add_argument(
        "--selection", choices=("all", "polymers-and-ligands", "polymers"), default="all"
    )
    parser.add_argument("--bonds", choices=("none", "explicit", "templates", "hybrid"))
    parser.add_argument("--conformers", choices=("first", "all"), default="all")
    parser.add_argument("--classifier", choices=CLASSIFIERS, help="Count atoms by classifier type")
    arguments = parser.parse_args()
    inspect_molecules(
        arguments.input,
        input_format=arguments.format,
        selection=arguments.selection,
        bonds=arguments.bonds,
        conformers=arguments.conformers,
        classifier=arguments.classifier,
    )


if __name__ == "__main__":
    main()
