"""Executable documentation example checks."""

from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path

import gemmi
import pytest

PROJECT_ROOT = Path(__file__).parents[2]
EXAMPLES = PROJECT_ROOT / "docs" / "examples"

WATER_PDB = """\
HETATM    1  O   HOH A   1       0.000   0.000   0.000  1.00 20.00           O
HETATM    2  H1  HOH A   1       0.957   0.000   0.000  1.00 20.00           H
HETATM    3  H2  HOH A   1      -0.240   0.927   0.000  1.00 20.00           H
CONECT    1    2    3
END
"""

WATER_MOL = """water
  ChargeFW          3D

  3  2  0  0  0  0  0  0  0  0  1 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.9570    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2400    0.9270    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
M  END
"""


WATERS_SDF = f"{WATER_MOL}$$$$\n{WATER_MOL}$$$$\n"


def run_example(name: str, *arguments: str | Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(EXAMPLES / name), *map(str, arguments)],
        check=True,
        capture_output=True,
        text=True,
    )


def test_calculate_file_example(tmp_path: Path) -> None:
    input_path = tmp_path / "waters.sdf"
    output_path = tmp_path / "result.json"
    csv_path = tmp_path / "charges.csv"
    input_path.write_text(WATERS_SDF, encoding="utf-8")

    run_example(
        "calculate_file.py",
        input_path,
        "--format",
        "sdf",
        "--method",
        "eem",
        "--result-json",
        output_path,
        "--csv",
        csv_path,
    )

    result = json.loads(output_path.read_text(encoding="utf-8"))
    assert result["status"] == "success"
    assert len(result["results"]) == 2
    with csv_path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    assert [(row["molecule"], row["element"]) for row in rows] == [
        (molecule, element) for molecule in ("0", "1") for element in ("O", "H", "H")
    ]
    assert all(row["conformer"] == "0" for row in rows)


def test_dipole_moment_example(tmp_path: Path) -> None:
    input_path = tmp_path / "waters.sdf"
    input_path.write_text(WATERS_SDF, encoding="utf-8")

    run_example("calculate_dipole_moments.py", input_path, "--format", "sdf", "--method", "eem")


def test_inspect_molecules_example(tmp_path: Path) -> None:
    input_path = tmp_path / "water.sdf"
    input_path.write_text(f"{WATER_MOL}$$$$\n", encoding="utf-8")

    run_example("inspect_molecules.py", input_path, "--format", "sdf", "--classifier", "hbo")


def test_gemmi_document_example_preserves_unchanged_mmcif(tmp_path: Path) -> None:
    input_path = tmp_path / "water.cif"
    output_path = tmp_path / "charged.cif"
    document = gemmi.read_pdb_string(WATER_PDB).make_mmcif_document()
    document.sole_block().set_pair("_audit.creation_method", "recipe-test")
    document.write_file(str(input_path))

    run_example("charge_gemmi_document.py", input_path, output_path, "--format", "mmcif")

    block = gemmi.cif.read_file(str(output_path)).sole_block()
    assert block.find_value("_audit.creation_method") == "recipe-test"
    assert "_sb_ncbr_partial_atomic_charges." in block.get_mmcif_category_names()


def test_fixed_ion_example(tmp_path: Path) -> None:
    input_path = PROJECT_ROOT / "tests" / "fixtures" / "synthetic" / "cif" / "fixed_ions.cif"
    output_path = tmp_path / "charges.json"

    run_example(
        "calculate_with_fixed_ions.py",
        input_path,
        "--format",
        "mmcif",
        "--ion",
        "MG",
        "--ion",
        "CA",
        "--result-json",
        output_path,
    )

    result = json.loads(output_path.read_text(encoding="utf-8"))
    assert result["status"] == "success"
    assert result["calculation_provenance"]["effective"]["fixed_ions"][0]["component_id"] == "MG"


def test_rdkit_sdf_example(tmp_path: Path) -> None:
    chemistry = pytest.importorskip("rdkit.Chem")
    input_path = tmp_path / "waters.sdf"
    output_path = tmp_path / "charged.sdf"
    input_path.write_text(WATERS_SDF, encoding="utf-8")

    run_example("charge_rdkit_sdf.py", input_path, output_path, "--method", "eem")

    molecules = list(chemistry.SDMolSupplier(str(output_path), removeHs=False))
    assert len(molecules) == 2
    assert all(
        atom.HasProp("ChargeFWPartialCharge")
        for molecule in molecules
        for atom in molecule.GetAtoms()
    )


def test_rdkit_conformer_example() -> None:
    pytest.importorskip("rdkit")
    run_example("analyze_rdkit_conformer_charges.py", "CCO", "--conformers", "3")


def test_compare_parameter_sets_example(tmp_path: Path) -> None:
    input_path = tmp_path / "water.sdf"
    input_path.write_text(f"{WATER_MOL}$$$$\n", encoding="utf-8")

    run_example("compare_parameter_sets.py", input_path, "--format", "sdf", "--method", "qeq")
