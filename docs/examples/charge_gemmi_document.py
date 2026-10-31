"""Attach calculated charges to a Gemmi mmCIF document.

Run with: python charge_gemmi_document.py input.cif charged.cif
"""

import argparse
from pathlib import Path

import chargefw
import chargefw.io.gemmi
import gemmi


def charge_document(
    input_path: str | Path,
    output_path: str | Path,
    *,
    method: str = "eem",
    parameter_set: str | None = None,
    selection: str = "all",
    result_json: str | Path | None = None,
) -> None:
    """Calculate charges and add them to the mmCIF document, keeping all its other data."""

    input_path = Path(input_path)
    document = gemmi.cif.read_file(str(input_path))
    molecules = chargefw.io.gemmi.from_document(
        document,
        source_name=str(input_path),
        selection=selection,
    )
    result = chargefw.calculate(
        molecules,
        method=method,
        parameter_set=parameter_set,
    )
    for warning in result.warnings:
        print(f"Warning: {warning.message}")

    chargefw.io.gemmi.attach_charges(document, result)
    document.write_file(str(output_path))
    if result_json is not None:
        chargefw.io.write(result_json, result, format="result-json")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--method", default="eem")
    parser.add_argument("--parameter-set")
    parser.add_argument(
        "--selection", choices=("all", "polymers-and-ligands", "polymers"), default="all"
    )
    parser.add_argument("--result-json", type=Path)
    arguments = parser.parse_args()
    charge_document(
        arguments.input,
        arguments.output,
        method=arguments.method,
        parameter_set=arguments.parameter_set,
        selection=arguments.selection,
        result_json=arguments.result_json,
    )


if __name__ == "__main__":
    main()
