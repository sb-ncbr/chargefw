"""Generated molecular and calculation-result output checks."""

import json
from gc import collect
from pathlib import Path
from typing import Any, cast

import chargefw
import gemmi
import jsonschema
import numpy as np
import pytest

SCHEMAS = Path(__file__).parents[2] / "schemas"
AROMATIC_MOL2 = Path(__file__).parents[1] / "fixtures" / "synthetic" / "mol2" / "aromatic.mol2"
TOPOLOGY_JSON = """{
  "schema_version": "1.0",
  "molecules": [{
    "atoms": [
      {"atomic_number": 8, "formal_charge": 0},
      {"atomic_number": 1, "formal_charge": 0},
      {"atomic_number": 1, "formal_charge": 0}
    ],
    "bonds": [
      {"atoms": [0, 1], "order": 1},
      {"atoms": [0, 2], "order": 1}
    ]
  }]
}"""


def validate(instance: Any, schema_name: str) -> None:
    schema = json.loads((SCHEMAS / schema_name).read_text(encoding="utf-8"))
    jsonschema.validate(instance, schema)


def result_json(result: chargefw.CalculationResult) -> dict[str, Any]:
    encoded: dict[str, Any] = json.loads(chargefw.io.dumps(result, format="result-json"))
    validate(encoded, "result-1.0.schema.json")
    return encoded


def water(*, conformers: int = 1) -> chargefw.Molecule:
    coordinates = [
        [[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]],
        [[0.1, 0.0, 0.0], [1.06, 0.0, 0.0], [-0.14, 0.93, 0.0]],
    ][:conformers]
    return chargefw.Molecule(
        [8, 1, 1],
        bonds=[[0, 1, 1], [0, 2, 1]],
        coordinates=coordinates,
        name="water",
        source_name="water-input",
        record_id="water-1",
    )


def test_formats_use_native_writers() -> None:
    molecule = water()
    result = chargefw.calculate(molecule, method="formal")

    mol2 = chargefw.io.dumps(result, format="mol2")
    assert "@<TRIPOS>MOLECULE\nwater-1" in mol2
    assert "USER_CHARGES" in mol2

    mmcif = chargefw.io.dumps(result, format="mmcif")
    assert "data_water-1" in mmcif
    assert "_sb_ncbr_partial_atomic_charges." in mmcif

    encoded = result_json(result)
    assert encoded["status"] == "success"
    assert encoded["results"][0]["input"]["source"] == "water-input"
    assert encoded["results"][0]["input"]["record_id"] == "water-1"
    assert encoded["calculation_provenance"]["effective"]["method"]["id"] == "formal"
    assert result.molecules[0] is molecule


def test_pqr_requires_structural_input() -> None:
    structure = chargefw.io.parse(
        "ATOM      1  O   HOH A   1       0.000   0.000   0.000  1.00 20.00           O  \n"
        "ATOM      2  H1  HOH A   1       0.957   0.000   0.000  1.00 20.00           H  \n"
        "ATOM      3  H2  HOH A   1      -0.239   0.927   0.000  1.00 20.00           H  \n"
        "END\n",
        format="pdb",
    )
    pqr = chargefw.io.dumps(chargefw.calculate(structure, method="formal"), format="pqr")
    assert "ATOM      2  H1  HOH A    1" in pqr

    with pytest.raises(ValueError, match="PDB or mmCIF"):
        chargefw.io.dumps(chargefw.calculate(water(), method="formal"), format="pqr")


def test_mmcif_holds_one_charge_set_per_result() -> None:
    molecule = water(conformers=2)
    formal = chargefw.calculate(molecule, method="formal")
    qeq = chargefw.calculate(molecule, method="qeq")
    block = gemmi.cif.read_string(chargefw.io.dumps([formal, qeq], format="mmcif")).sole_block()
    methods = block.find("_sb_ncbr_partial_atomic_charges_meta.", ["id", "method"])
    assert [(gemmi.cif.as_string(row[0]), gemmi.cif.as_string(row[1])) for row in methods] == [
        ("1", "formal"),
        ("2", "qeq"),
    ]
    # Each set covers both conformers under its single id.
    type_ids = [
        gemmi.cif.as_string(row[0])
        for row in block.find("_sb_ncbr_partial_atomic_charges.", ["type_id"])
    ]
    assert type_ids == ["1"] * 6 + ["2"] * 6

    with pytest.raises(ValueError, match="one calculation result"):
        chargefw.io.dumps([formal, qeq], format="mol2")
    for rejected in (
        [formal, chargefw.calculate(water(), method="qeq")],
        [formal, chargefw.calculate(molecule, method="formal")],
    ):
        with pytest.raises(ValueError):
            chargefw.io.dumps(rejected, format="mmcif")


def test_write_uses_explicit_format(tmp_path: Path) -> None:
    result = chargefw.calculate(water(), method="formal")
    path = tmp_path / "charges.data"
    chargefw.io.write(path, result, format="mol2")
    assert "@<TRIPOS>MOLECULE" in path.read_text(encoding="utf-8")


@pytest.mark.parametrize(("selection", "expected_count"), [("all", 2), ("first", 1)])
def test_mmcif_applies_geometry_independent_charges_to_selected_conformers(
    selection: str, expected_count: int
) -> None:
    contents = json.dumps(
        {
            "schema_version": "1.0",
            "molecules": [
                {
                    "atoms": [{"atomic_number": 1, "formal_charge": 0}],
                    "conformers": [
                        {"coordinates": [[0, 0, 0]]},
                        {"coordinates": [[1, 0, 0]]},
                    ],
                }
            ],
        }
    )
    molecules = chargefw.io.parse(contents, format="molecule-json", conformers=cast(Any, selection))
    result = chargefw.calculate(molecules, method="formal")
    block = gemmi.cif.read_string(chargefw.io.dumps(result, format="mmcif")).sole_block()
    atom_sites = block.find("_atom_site.", ["id"])
    charges = block.find("_sb_ncbr_partial_atomic_charges.", ["atom_id"])
    assert len(atom_sites) == expected_count
    assert [gemmi.cif.as_string(row[0]) for row in charges] == [
        gemmi.cif.as_string(row[0]) for row in atom_sites
    ]


def test_failed_result_json_preserves_input_identity() -> None:
    validate(json.loads(TOPOLOGY_JSON), "molecule-input-1.0.schema.json")
    imported = chargefw.io.parse(TOPOLOGY_JSON, format="molecule-json", source_name="topology.json")
    manual = chargefw.Molecule([1], atom_ids=["hydrogen"])
    with pytest.raises(chargefw.NoExecutablePlanError) as raised:
        chargefw.calculate([imported[0], manual], method="qeq")
    failed_result = raised.value.result

    encoded = result_json(failed_result)
    assert encoded["status"] == "no_executable_plan"
    imported_record, manual_record = encoded["results"]
    assert imported_record["input"]["import"]["format"] == "molecule-json"
    assert manual_record["input"]["atom_ids"] == ["hydrogen"]
    assert all("assignments" not in record for record in encoded["results"])
    assert imported_record["diagnostics"]
    for format_name in ("mol2", "mmcif"):
        with pytest.raises(ValueError, match="successful calculation"):
            chargefw.io.dumps(failed_result, format=cast(Any, format_name))


def test_import_metadata_follows_reordered_and_recombined_molecules() -> None:
    imported = chargefw.io.read(AROMATIC_MOL2, format="mol2")
    imported_molecule = imported[0]
    del imported
    collect()

    manual = chargefw.Molecule([1], source_name="manual")
    result = chargefw.calculate([manual, imported_molecule], method="formal")
    manual_record, imported_record = result_json(result)["results"]

    assert manual_record["input"] == {"source": "manual", "record_index": 0}
    assert manual_record["diagnostics"] == []
    assert imported_record["input"]["import"]["format"] == "mol2"
    assert "atom_ids" not in imported_record["input"]
    assert len(imported_record["assignments"][0]["charges"]) == len(imported_molecule.atom_ids)
    assert [value["code"] for value in imported_record["diagnostics"]] == [
        "partial_charges_ignored"
    ]


def test_coordinate_free_result_json_is_molecule_scoped() -> None:
    molecules = chargefw.io.parse(
        """{
  "schema_version": "1.0",
  "molecules": [{"atoms": [{"atomic_number": 1, "formal_charge": 0}]}]
}""",
        format="molecule-json",
        source_name="hydrogen.json",
    )
    result = chargefw.calculate(molecules, method="formal")
    del molecules
    collect()

    record = result_json(result)["results"][0]
    assert record["assignments"][0]["scope"] == "molecule"
    assert record["input"]["import"]["format"] == "molecule-json"


def test_result_outputs_preserve_integer_caller_ids() -> None:
    molecule = chargefw.Molecule(
        [1, 1],
        coordinates=[[0, 0, 0], [1, 0, 0]],
        name="hydrogen",
        record_id=42,
        atom_ids=cast(Any, ["left", np.int64(9)]),
    )
    result = chargefw.calculate(molecule, method="formal")

    encoded = result_json(result)
    assert encoded["results"][0]["input"]["record_id"] == 42
    assert encoded["results"][0]["input"]["atom_ids"] == ["left", 9]
    assert result.assignments[0].atom_ids == ("left", 9)
    assert "\n42\n" in chargefw.io.dumps(result, format="mol2")
    assert "data_42" in chargefw.io.dumps(result, format="mmcif")


@pytest.mark.parametrize("format_name", ["mol2", "mmcif"])
def test_molecular_output_requires_finite_coordinates(format_name: Any) -> None:
    missing = chargefw.calculate(chargefw.Molecule([1]), method="formal")
    result_json_text = chargefw.io.dumps(missing, format="result-json")
    with pytest.raises(ValueError, match="coordinates"):
        chargefw.io.dumps(missing, format=format_name)
    assert chargefw.io.dumps(missing, format="result-json") == result_json_text

    nonfinite = chargefw.calculate(
        chargefw.Molecule([1], coordinates=[[float("nan"), 0.0, 0.0]]),
        method="formal",
    )
    with pytest.raises(ValueError, match="finite"):
        chargefw.io.dumps(nonfinite, format=format_name)


def test_output_format_is_explicit() -> None:
    result = chargefw.calculate(water(), method="formal")
    with pytest.raises(TypeError):
        chargefw.io.dumps(result)  # type: ignore[call-arg]
