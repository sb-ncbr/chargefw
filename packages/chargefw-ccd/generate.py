#!/usr/bin/env python3
"""Writes the chargefw-ccd data file from the wwPDB Chemical Component Dictionary.

Usage: generate.py components.cif.gz

Writes two files into src/chargefw_ccd/:

    components.txt.gz   the atoms and bonds of every component
    _version.py         the package version: the date of the newest component, followed by
                        .postREVISION when REVISION below is not zero

The data file is gzip-compressed text. Its first line is

    chargefw-ccd FORMAT_VERSION PACKAGE_VERSION

Each component follows as a line

    ID ATOM_COUNT BOND_COUNT

then one line per atom,

    NAME ELEMENT LEAVING

where LEAVING is 1 for an atom that is absent when the component is linked in a polymer and 0
otherwise, then one line per bond,

    FIRST SECOND ORDER

where FIRST and SECOND are zero-based positions in the component's atom list and ORDER is 1, 2, or 3.
The same dictionary always gives the same two files, byte for byte.
"""

import gzip
import sys
from pathlib import Path

import gemmi
from gemmi import cif

FORMAT_VERSION = 1
# A dictionary date can be published only once. To publish a corrected file for a date that is already
# released, raise this number.
REVISION = 0
BOND_ORDERS = {"SING": 1, "DOUB": 2, "TRIP": 3}
LEAVING = {"N": 0, "Y": 1}
PACKAGE = Path(__file__).parent / "src" / "chargefw_ccd"


def rows(block: cif.Block, category: str, columns: list[str]) -> cif.Table:
    """Returns the given columns of a category, which the component may lack as a whole."""
    table = block.find(category, columns)
    if table.width() == 0 and block.find_mmcif_category(category).width() != 0:
        raise ValueError(f"{block.name}: {category} lacks one of {', '.join(columns)}")
    return table


def component_lines(block: cif.Block) -> list[str]:
    """Returns the lines of one component."""
    name = block.name
    if len(name.split()) != 1:
        raise ValueError(f"component name {name!r} is not a single word")

    atoms: dict[str, int] = {}
    lines = []
    for row in rows(
        block, "_chem_comp_atom.", ["atom_id", "type_symbol", "pdbx_leaving_atom_flag"]
    ):
        atom, symbol, leaving = (row.str(column) for column in range(3))
        if len(atom.split()) != 1:
            raise ValueError(f"{name}: atom name {atom!r} is not a single word")
        if atom in atoms:
            raise ValueError(f"{name}: atom {atom} is listed twice")
        element = gemmi.Element(symbol)
        if element.atomic_number == 0 and symbol != "X":
            raise ValueError(f"{name}: atom {atom} has unknown element {symbol!r}")
        if leaving not in LEAVING:
            raise ValueError(f"{name}: atom {atom} has leaving flag {leaving!r}")
        atoms[atom] = len(atoms)
        lines.append(f"{atom} {element.name} {LEAVING[leaving]}")

    bond_count = 0
    for row in rows(block, "_chem_comp_bond.", ["atom_id_1", "atom_id_2", "value_order"]):
        first, second, order = (row.str(column) for column in range(3))
        for atom in (first, second):
            if atom not in atoms:
                raise ValueError(f"{name}: bond to atom {atom}, which the component does not have")
        if order not in BOND_ORDERS:
            raise ValueError(f"{name}: bond {first}-{second} has order {order!r}")
        lines.append(f"{atoms[first]} {atoms[second]} {BOND_ORDERS[order]}")
        bond_count += 1

    return [f"{name} {len(atoms)} {bond_count}", *lines]


def main() -> None:
    if len(sys.argv) != 2:
        sys.exit(__doc__.split("\n\n")[1])
    document = cif.read(sys.argv[1])

    lines = []
    dates = []
    for block in document:
        lines += component_lines(block)
        dates.append(cif.as_string(block.find_value("_chem_comp.pdbx_modified_date")))
    year, month, day = (int(part) for part in max(dates).split("-"))

    version = f"{year}.{month}.{day}" + (f".post{REVISION}" if REVISION else "")
    text = "\n".join([f"chargefw-ccd {FORMAT_VERSION} {version}", *lines, ""])
    # No file name or time in the gzip header, so that the output depends only on the dictionary.
    with (
        open(PACKAGE / "components.txt.gz", "wb") as file,
        gzip.GzipFile(filename="", mode="wb", fileobj=file, mtime=0) as compressed,
    ):
        compressed.write(text.encode("ascii"))
    (PACKAGE / "_version.py").write_text(f'__version__ = "{version}"\n')
    print(version)


if __name__ == "__main__":
    main()
