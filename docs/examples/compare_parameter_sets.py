"""Assess and compare bundled parameter sets for one method and one molecular geometry.

Run with: python compare_parameter_sets.py molecule.sdf --format sdf --method eem
"""

import argparse
from pathlib import Path

import chargefw
import numpy as np


def compare_parameter_sets(
    input_path: str | Path,
    *,
    input_format: str,
    method_id: str,
    reference_id: str | None = None,
) -> None:
    """Assess once and report differences without ranking scientific accuracy."""

    molecules = chargefw.io.read(input_path, format=input_format)
    if len(molecules) != 1 or molecules[0].conformer_count != 1:
        raise ValueError("this example requires exactly one molecule with one conformer")

    method = chargefw.methods[method_id]
    if not method.parameter_sets:
        raise ValueError(f"method {method_id!r} has no bundled parameter sets")

    assessment = chargefw.assess(molecules, method=method, execution="full")
    for rejection in assessment.rejections:
        parameter_set = rejection.parameter_set
        identifier = parameter_set.id if parameter_set is not None else "-"
        reasons = "; ".join(issue.message for issue in rejection.issues)
        print(f"Rejected {identifier}: {reasons}")

    executable: list[tuple[str, np.ndarray]] = []
    for plan in assessment.plans:
        parameter_set = plan.parameter_set
        if parameter_set is None:
            raise RuntimeError(f"{method_id!r} produced a plan without a parameter set")
        result = chargefw.calculate(molecules, plan)
        executable.append((parameter_set.id, result.assignments[0].values))

    if not executable:
        raise RuntimeError(f"no bundled {method_id} parameter set is applicable")

    if reference_id is None:
        reference_id = executable[0][0]
    try:
        reference = next(values for identifier, values in executable if identifier == reference_id)
    except StopIteration as error:
        raise ValueError(f"reference parameter set {reference_id!r} was not executable") from error

    print(f"Reference: {reference_id}")
    print(
        "Parameter set                              total charge    RMS difference    max difference"
    )
    for identifier, values in executable:
        difference = values - reference
        rms = float(np.sqrt(np.mean(np.square(difference))))
        maximum = float(np.max(np.abs(difference)))
        print(f"{identifier:<42} {values.sum():12.6f} {rms:17.6f} {maximum:17.6f}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--format", choices=chargefw.io.INPUT_FORMATS, required=True)
    parser.add_argument("--method", default="eem")
    parser.add_argument("--reference")
    arguments = parser.parse_args()
    compare_parameter_sets(
        arguments.input,
        input_format=arguments.format,
        method_id=arguments.method,
        reference_id=arguments.reference,
    )


if __name__ == "__main__":
    main()
