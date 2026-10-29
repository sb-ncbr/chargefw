"""Molecular format input and Gemmi conversion checks."""

import json
from collections.abc import Sequence
from dataclasses import FrozenInstanceError
from pathlib import Path
from typing import Any, cast, get_type_hints

import chargefw.io.gemmi
import gemmi
import numpy as np
import pytest
from chargefw import calculate
from chargefw import io as chargefw_io

RESULTS = chargefw.CalculationResult | Sequence[chargefw.CalculationResult]
PDB_TEXT = """HEADER    TEST PDB
TITLE     TWO MODELS
MODEL        1
ATOM      1  O   HOH A   1       0.000   0.000   0.000  1.00 20.00           O  
ATOM      2  H1 AHOH A   1       0.957   0.000   0.000  1.00 20.00           H  
ATOM      3  H1 BHOH A   1       9.000   0.000   0.000  1.00 20.00           H  
ENDMDL
MODEL        2
ATOM      1  O   HOH A   1       0.100   0.000   0.000  1.00 20.00           O  
ATOM      2  H1 AHOH A   1       1.057   0.000   0.000  1.00 20.00           H  
ATOM      3  H1 BHOH A   1       9.100   0.000   0.000  1.00 20.00           H  
ENDMDL
CONECT    1    2
END
"""

MMCIF_TEXT = """data_first
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 C CA . ALA A 1 ? 0.0 0.0 0.0 1.0 20.0 0 1 ALA A CA 1
HETATM 2 O O . HOH A 2 ? 1.0 0.0 0.0 1.0 20.0 0 2 HOH A O 1
ATOM 3 C CA . ALA A 1 ? 0.1 0.0 0.0 1.0 20.0 0 1 ALA A CA 2
HETATM 4 O O . HOH A 2 ? 1.1 0.0 0.0 1.0 20.0 0 2 HOH A O 2
#
data_second
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 O O . HOH B 1 ? 3.0 0.0 0.0 1.0 20.0 0 1 HOH B O 1
#
"""

FIXED_IONS_MMCIF_TEXT = """data_fixed_ions
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.label_entity_id
_atom_site.pdbx_PDB_model_num
HETATM 1 O O . HOH W . ? 0 0 0 1 20 0 5 HOH W O 1 1
HETATM 2 H H1 . HOH W . ? 0.96 0 0 1 20 0 5 HOH W H1 1 1
HETATM 3 H H2 . HOH W . ? -0.24 0.93 0 1 20 0 5 HOH W H2 1 1
HETATM 4 Mg MG . MG I . ? 0 0 3 1 20 0 9 AUTH_MG I AUTH_MG 2 1
HETATM 5 O O . HOH W . ? 0 0 0.1 1 20 0 5 HOH W O 1 2
HETATM 6 H H1 . HOH W . ? 0.96 0 0.1 1 20 0 5 HOH W H1 1 2
HETATM 7 H H2 . HOH W . ? -0.24 0.93 0.1 1 20 0 5 HOH W H2 1 2
HETATM 8 Mg MG . MG I . ? 0.25 0 3 1 20 0 9 AUTH_MG I AUTH_MG 2 2
#
"""

BOND_STRATEGY_PDB = """SSBOND   1 CYS A   3    CYS A   4                          
LINK         C   ALA A   1                 C1  LIG A   5
ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N  
ATOM      2  CA  ALA A   1       1.450   0.000   0.000  1.00 20.00           C  
ATOM      3  C   ALA A   1       2.450   1.000   0.000  1.00 20.00           C  
ATOM      4  O   ALA A   1       3.450   1.000   0.000  1.00 20.00           O  
ATOM      5  CB  ALA A   1       1.450  -1.000   0.000  1.00 20.00           C  
ATOM      6  N   GLY A   2       2.200   2.200   0.000  1.00 20.00           N  
ATOM      7  CA  GLY A   2       3.200   3.200   0.000  1.00 20.00           C  
ATOM      8  C   GLY A   2       4.200   3.200   0.000  1.00 20.00           C  
ATOM      9  O   GLY A   2       5.200   3.200   0.000  1.00 20.00           O  
ATOM     10  SG  CYS A   3       6.200   3.200   0.000  1.00 20.00           S  
ATOM     11  SG  CYS A   4       7.200   3.200   0.000  1.00 20.00           S  
HETATM   12  C1  LIG A   5       8.200   3.200   0.000  1.00 20.00           C  
CONECT   11   12
END
"""

MOL_TEXT = """charged-v2000
  ChargeFW

  2  1  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
    1.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
M  CHG  2   1   1   2  -1
M  END
"""

MOL2_TEXT = """@<TRIPOS>MOLECULE
water
3 2 0 0 0
SMALL
NO_CHARGES

@<TRIPOS>ATOM
1 O 0.0000 0.0000 0.0000 O.2
2 H1 0.9570 0.0000 0.0000 H
3 H2 -0.2400 0.9270 0.0000 H
@<TRIPOS>BOND
1 1 2 1
2 1 3 1
"""

MOLECULE_JSON_TEXT = """{
  "schema_version": "1.0",
  "molecules": [{
    "id": "water-1",
    "name": "water",
    "atoms": [
      {"atomic_number": 8, "formal_charge": 0},
      {"atomic_number": 1, "formal_charge": 0}
    ],
    "bonds": [{"atoms": [0, 1], "order": 1}],
    "conformers": [
      {"id": "first", "coordinates": [[0, 0, 0], [1, 0, 0]]},
      {"id": "second", "coordinates": [[0.1, 0, 0], [1.1, 0, 0]]}
    ]
  }]
}
"""


def test_fixed_ions_resolve_gemmi_components_and_provenance() -> None:
    molecules = chargefw_io.parse(
        FIXED_IONS_MMCIF_TEXT,
        format="mmcif",
        conformers=cast(Any, "all"),
        bonds="templates",
    )
    assert molecules[0].formal_charges.tolist() == [0, 0, 0, 0]
    mapping = molecules[0].source_mapping
    assert mapping is not None
    labels = mapping.atoms[3].structural_labels
    assert labels is not None
    assert labels.label.residue == "MG"
    assert labels.author.residue == "AUTH_MG"
    result = chargefw.calculate(
        molecules,
        method="sqeqp",
        parameter_set="SQEqp_Schindler2021_CCD_gen",
        execution="full",
        fixed_ions=["MG"],
    )

    encoded = json.loads(chargefw_io.dumps(result, format="result-json"))
    effective = encoded["calculation_provenance"]["effective"]["fixed_ions"]
    assert effective == [
        {
            "component_id": "MG",
            "charge": 2.0,
            "instances": [{"molecule_index": 0, "atom_index": 3}],
        }
    ]


def test_parse_native_molecular_formats() -> None:
    molecules = chargefw_io.parse(MOL_TEXT, format="mol", source_name="charged.mol")
    assert len(molecules) == 1
    molecule = molecules[0]
    assert molecule.atomic_numbers.tolist() == [7, 8]
    assert molecule.formal_charges.tolist() == [1, -1]
    assert molecule.source_name == "charged.mol"

    sdf = chargefw_io.parse(f"{MOL_TEXT}$$$$\n", format="sdf", source_name="charged.sdf")
    assert len(sdf) == 1
    assert sdf[0].record_index == 0

    mol2 = chargefw_io.parse(MOL2_TEXT, format="mol2", source_name="water.mol2")
    assert len(mol2) == 1
    assert mol2[0].atom_names == ("O", "H1", "H2")
    assert mol2[0].bond_count == 2


def test_parse_molecule_json_conformer_selection() -> None:
    all_conformers = chargefw_io.parse(MOLECULE_JSON_TEXT, format="molecule-json")
    first_conformer = chargefw_io.parse(
        MOLECULE_JSON_TEXT, format="molecule-json", conformers="first"
    )

    assert all_conformers[0].conformer_names == ("first", "second")
    assert first_conformer[0].conformer_names == ("first",)


def test_parse_topology_only_molecule_json_supports_formal_calculation() -> None:
    topology_only_conformers: tuple[list[Any] | None, ...] = (None, [])
    for conformers in topology_only_conformers:
        molecule_data: dict[str, Any] = {
            "schema_version": "1.0",
            "molecules": [
                {
                    "atoms": [
                        {"atomic_number": 7, "formal_charge": 1},
                        {"atomic_number": 1, "formal_charge": 0},
                    ],
                    "bonds": [{"atoms": [0, 1], "order": 1}],
                }
            ],
        }
        if conformers is not None:
            molecule_data["molecules"][0]["conformers"] = conformers

        molecules = chargefw_io.parse(json.dumps(molecule_data), format="molecule-json")
        assert molecules[0].conformer_count == 0
        result = calculate(molecules, method="formal", execution="full")

        assert result.status == "success"
        np.testing.assert_array_equal(result.assignments[0].values, [1.0, 0.0])


def test_path_readers_set_source_name(tmp_path: Path) -> None:
    path = tmp_path / "charged.mol"
    path.write_text(MOL_TEXT, encoding="utf-8")
    molecules = chargefw_io.read(path, format="mol")

    assert molecules[0].source_name == str(path)


def test_public_annotations_resolve_to_runtime_types() -> None:
    adapter = chargefw.io.gemmi
    cases = (
        (
            adapter.from_structure,
            {"structure": gemmi.Structure, "return": chargefw.MoleculeCollection},
        ),
        (
            adapter.from_document,
            {"document": gemmi.cif.Document, "return": chargefw.MoleculeCollection},
        ),
        (
            adapter.to_document,
            {"result": RESULTS, "return": gemmi.cif.Document},
        ),
        (
            adapter.attach_charges,
            {
                "document": gemmi.cif.Document,
                "result": RESULTS,
                "return": type(None),
            },
        ),
    )
    for function, expected in cases:
        hints = get_type_hints(function)
        for name, expected_type in expected.items():
            assert hints[name] == expected_type


def test_object_conversion_reports_missing_optional_dependency(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def missing(name: str) -> Any:
        raise ModuleNotFoundError(f"No module named {name!r}", name=name)

    monkeypatch.setattr(chargefw.io.gemmi, "import_module", missing)
    with pytest.raises(ImportError, match=r"pip install chargefw\[gemmi\]"):
        chargefw.io.gemmi.from_structure(cast(Any, object()))


def test_pdb_text_file_and_structure_preserve_mapping(tmp_path: Path) -> None:
    text_collection = chargefw_io.parse(
        PDB_TEXT,
        format="pdb",
        source_name="water.pdb",
        bonds="explicit",
    )
    assert len(text_collection) == 1
    text_molecule = text_collection[0]
    assert text_molecule.atom_names == ("O", "H1")
    assert text_molecule.atom_ids == ("1", "2")
    assert text_molecule.conformer_names == ("1", "2")
    assert text_molecule.source_name == "water.pdb"
    assert text_molecule.bonds.tolist() == [[0, 1, 1]]
    np.testing.assert_allclose(text_molecule.coordinates[:, 0, 0], [0.0, 0.1])
    mapping = text_molecule.source_mapping
    assert mapping is not None
    assert mapping.format == "pdb"
    assert mapping.components == (chargefw.SourceComponentInstance("HOH", (0, 1)),)
    assert mapping.alternate_location_selection == "first-source-order"
    assert [value.id for value in mapping.conformers] == ["1", "2"]
    labels = mapping.atoms[1].structural_labels
    assert labels is not None
    assert labels.author.atom == "H1"
    assert labels.author.residue == "HOH"
    assert labels.author.chain == "A"
    assert labels.alternate_location == "A"

    path = tmp_path / "water.pdb"
    path.write_text(PDB_TEXT, encoding="utf-8")
    file_molecule = chargefw_io.read(path, format="pdb", bonds="explicit")[0]
    np.testing.assert_array_equal(file_molecule.atomic_numbers, text_molecule.atomic_numbers)
    np.testing.assert_array_equal(file_molecule.coordinates, text_molecule.coordinates)

    structure = gemmi.read_pdb_string(PDB_TEXT)
    structure_collection = chargefw.io.gemmi.from_structure(structure, source_name="structure")
    assert len(structure_collection) == 1
    structure_molecule = structure_collection[0]
    np.testing.assert_array_equal(structure_molecule.atomic_numbers, text_molecule.atomic_numbers)
    np.testing.assert_array_equal(structure_molecule.coordinates, text_molecule.coordinates)


def test_mmcif_text_file_and_document_preserve_records(tmp_path: Path) -> None:
    collection = chargefw_io.parse(MMCIF_TEXT, format="mmcif", source_name="models.cif")
    assert len(collection) == 2
    assert [value.record_index for value in collection] == [0, 1]
    assert [value.record_id for value in collection] == ["first", "second"]
    assert collection[0].conformer_names == ("1", "2")
    assert collection[0].atom_ids == ("1", "2")
    mapping = collection[0].source_mapping
    assert mapping is not None
    assert [value.id for value in mapping.conformers] == ["1", "2"]
    assert [site.id for site in mapping.conformers[1].sites] == ["3", "4"]
    labels = mapping.atoms[0].structural_labels
    assert labels is not None
    assert labels.author.atom == "CA"
    assert labels.label.atom == "CA"

    path = tmp_path / "models.cif"
    path.write_text(MMCIF_TEXT, encoding="utf-8")
    from_file = chargefw_io.read(path, format="mmcif")
    assert [value.atom_count for value in from_file] == [2, 1]

    document = gemmi.cif.read_string(MMCIF_TEXT)
    from_document = chargefw.io.gemmi.from_document(document, source_name="document")
    assert [value.record_id for value in from_document] == ["first", "second"]
    np.testing.assert_array_equal(from_document[0].coordinates, collection[0].coordinates)


def test_mmcif_mapping_exposes_arbitrary_ids_and_label_namespaces() -> None:
    contents = """data_mapping
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.label_entity_id
_atom_site.pdbx_PDB_model_num
HETATM 001 C LABEL . LIG LC 7 ? 0 0 0 1 20 0 17 AUTH AC AUTHOR E1 1
#
"""
    molecule = chargefw_io.parse(contents, format="mmcif")[0]
    del contents

    assert molecule.atom_ids == ("001",)
    mapping = molecule.source_mapping
    assert mapping is not None
    labels = mapping.atoms[0].structural_labels
    assert labels is not None
    assert labels.author.atom == "AUTHOR"
    assert labels.author.residue == "AUTH"
    assert labels.author.chain == "AC"
    assert labels.label.atom == "LABEL"
    assert labels.label.residue == "LIG"
    assert labels.label.chain == "LC"
    assert mapping.components == (chargefw.SourceComponentInstance("LIG", (0,)),)
    assert isinstance(mapping.components, tuple)
    assert isinstance(mapping.components[0].atom_indices, tuple)
    with pytest.raises(FrozenInstanceError):
        setattr(mapping.components[0], "component_id", "AUTHOR")

    result = calculate(molecule, method="formal")
    assert result.molecules[0].source_mapping == mapping
    imported = json.loads(chargefw_io.dumps(result, format="result-json"))["results"][0]["input"][
        "import"
    ]
    assert imported["format"] == "mmcif"


def test_to_document_creates_fresh_mapped_mmcif() -> None:
    source = gemmi.cif.read_string(MMCIF_TEXT)
    source[0].set_pair("_audit.creation_method", "source-only")
    molecules = chargefw.io.gemmi.from_document(source)
    result = calculate(molecules, method="formal")

    generated = chargefw.io.gemmi.to_document(result)

    assert len(generated) == 2
    assert source[0].find_value("_audit.creation_method") == "source-only"
    assert "_sb_ncbr_partial_atomic_charges." not in source[0].get_mmcif_category_names()
    assert not generated[0].find_value("_audit.creation_method")
    sites = generated[0].find(
        "_atom_site.",
        ["id", "label_atom_id", "auth_atom_id", "label_alt_id", "pdbx_PDB_model_num"],
    )
    assert len(sites) == 4
    assert [gemmi.cif.as_string(row[0]) for row in sites] == ["1", "2", "3", "4"]
    assert gemmi.cif.as_string(sites[0][1]) == "CA"
    assert gemmi.cif.as_string(sites[0][2]) == "CA"
    assert gemmi.cif.as_string(sites[0][4]) == "1"
    assert gemmi.cif.as_string(sites[2][4]) == "2"
    charges = generated[0].find(
        "_sb_ncbr_partial_atomic_charges.", ["type_id", "atom_id", "charge"]
    )
    assert len(charges) == 4
    assert [gemmi.cif.as_string(row[1]) for row in charges] == ["1", "2", "3", "4"]
    with pytest.raises(TypeError):
        chargefw.io.gemmi.to_document(cast(Any, object()))


def test_attach_charges_requires_the_unchanged_import_document() -> None:
    document = gemmi.cif.read_string(MMCIF_TEXT)
    document[0].set_pair("_audit.creation_method", "attachment-test")
    molecules = chargefw.io.gemmi.from_document(document)
    result = calculate(molecules, method="formal")

    chargefw.io.gemmi.attach_charges(document, result)

    assert document[0].find_value("_audit.creation_method") == "attachment-test"
    first_ids = [
        gemmi.cif.as_string(row[0])
        for row in document[0].find("_sb_ncbr_partial_atomic_charges.", ["atom_id"])
    ]
    assert first_ids == ["1", "2", "3", "4"]
    before_rejected_overwrite = document.as_string()
    with pytest.raises(ValueError):
        chargefw.io.gemmi.attach_charges(document, result)
    assert document.as_string() == before_rejected_overwrite
    chargefw.io.gemmi.attach_charges(document, result, overwrite=True)

    conformer_document = gemmi.cif.read_string(MMCIF_TEXT)
    conformer_result = calculate(chargefw.io.gemmi.from_document(conformer_document), method="qeq")
    assert [assignment.conformer_index for assignment in conformer_result.assignments] == [0, 1, 0]
    chargefw.io.gemmi.attach_charges(conformer_document, conformer_result)
    first_charge_rows = conformer_document[0].find(
        "_sb_ncbr_partial_atomic_charges.", ["type_id", "atom_id"]
    )
    assert [
        (gemmi.cif.as_string(row[0]), gemmi.cif.as_string(row[1])) for row in first_charge_rows
    ] == [("1", "1"), ("1", "2"), ("1", "3"), ("1", "4")]

    sets_document = gemmi.cif.read_string(MMCIF_TEXT)
    sets_molecules = chargefw.io.gemmi.from_document(sets_document)
    chargefw.io.gemmi.attach_charges(
        sets_document,
        [calculate(sets_molecules, method="formal"), calculate(sets_molecules, method="qeq")],
    )
    methods = sets_document[0].find("_sb_ncbr_partial_atomic_charges_meta.", ["id", "method"])
    assert [(gemmi.cif.as_string(row[0]), gemmi.cif.as_string(row[1])) for row in methods] == [
        ("1", "formal"),
        ("2", "qeq"),
    ]
    assert len(sets_document[0].find("_sb_ncbr_partial_atomic_charges.", ["atom_id"])) == 8

    # A rejected second result leaves the document as it was.
    unchanged_document = gemmi.cif.read_string(MMCIF_TEXT)
    unchanged_molecules = chargefw.io.gemmi.from_document(unchanged_document)
    before_rejected_sets = unchanged_document.as_string()
    with pytest.raises(ValueError):
        chargefw.io.gemmi.attach_charges(
            unchanged_document,
            [
                calculate(unchanged_molecules, method="formal"),
                calculate(unchanged_molecules, method="formal"),
            ],
        )
    assert unchanged_document.as_string() == before_rejected_sets

    altloc_document = gemmi.cif.read_string(MMCIF_TEXT.replace(" C CA . ALA", " C CA B ALA"))
    altloc_result = calculate(chargefw.io.gemmi.from_document(altloc_document), method="formal")
    chargefw.io.gemmi.attach_charges(altloc_document, altloc_result)
    assert (
        gemmi.cif.as_string(altloc_document[0].find("_atom_site.", ["label_alt_id"])[0][0]) == "B"
    )
    assert len(altloc_document[0].find("_sb_ncbr_partial_atomic_charges.", ["atom_id"])) == 4

    for column, value in (
        ("id", "99"),
        ("Cartn_x", "9.0"),
        ("pdbx_formal_charge", "1"),
        ("type_symbol", "O"),
        ("auth_atom_id", "CB"),
    ):
        modified = gemmi.cif.read_string(MMCIF_TEXT)
        modified_molecules = chargefw.io.gemmi.from_document(modified)
        modified_result = calculate(modified_molecules, method="formal")
        modified[0].find("_atom_site.", [column])[0][0] = value
        before_mismatch = modified.as_string()
        with pytest.raises(ValueError):
            chargefw.io.gemmi.attach_charges(modified, modified_result)
        assert modified.as_string() == before_mismatch
    with pytest.raises(TypeError):
        chargefw.io.gemmi.attach_charges(document, result, overwrite=cast(Any, 1))


def test_attach_charges_ignores_atom_site_categories_without_ids() -> None:
    document = gemmi.cif.read_string(
        f"{MMCIF_TEXT}data_auxiliary\nloop_\n_atom_site.label_atom_id\nX\n#\n"
    )
    result = calculate(chargefw.io.gemmi.from_document(document), method="formal")

    chargefw.io.gemmi.attach_charges(document, result)

    assert len(document[0].find("_sb_ncbr_partial_atomic_charges.", ["atom_id"])) == 4
    assert len(document[1].find("_sb_ncbr_partial_atomic_charges.", ["atom_id"])) == 1
    assert "_sb_ncbr_partial_atomic_charges." not in document[2].get_mmcif_category_names()


def test_selection_conformers_and_types_are_explicit() -> None:
    with pytest.raises(RuntimeError):
        chargefw_io.parse(MMCIF_TEXT, format="mmcif", selection="polymers", conformers="first")

    polymers = chargefw_io.parse(
        BOND_STRATEGY_PDB,
        format="pdb",
        selection="polymers",
        conformers="first",
    )
    assert polymers[0].conformer_count == 1
    result = calculate(polymers, method="formal")
    encoded = json.loads(chargefw_io.dumps(result, format="result-json"))
    policy = encoded["results"][0]["input"]["import"]["policy"]
    assert policy["conformer_selection"] == "first"
    assert policy["record_selection"] == "polymers"
    assert policy["bond_strategy"] == "hybrid"

    with pytest.raises(TypeError):
        chargefw.io.gemmi.from_structure(cast(Any, object()))
    with pytest.raises(ValueError):
        chargefw_io.parse(PDB_TEXT, format="pdb", selection=cast(Any, "invalid"))
    with pytest.raises(ValueError):
        chargefw_io.parse(
            PDB_TEXT,
            format="pdb",
            selection=cast(Any, "polymers_and_ligands"),
        )
    with pytest.raises(TypeError):
        chargefw_io.parse(PDB_TEXT, format="pdb", selection=cast(Any, 1))


def test_bond_strategy_is_forwarded_to_native_adapter() -> None:
    molecule = chargefw_io.parse(BOND_STRATEGY_PDB, format="pdb", bonds="explicit")[0]
    assert molecule.bond_count == 2


def test_source_atom_ids_distinguish_repeated_atom_names() -> None:
    molecule = chargefw_io.parse(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  \n"
        "ATOM      2  CA  GLY A   2       1.000   0.000   0.000  1.00 20.00           C  \n"
        "END\n",
        format="pdb",
    )[0]
    assert molecule.atom_names == ("CA", "CA")
    assert molecule.atom_ids == ("1", "2")


def test_generic_input_requires_explicit_compatible_format_options() -> None:
    with pytest.raises(TypeError):
        chargefw_io.parse(MOL_TEXT)  # type: ignore[call-arg]
    with pytest.raises(ValueError):
        chargefw_io.parse(MOL_TEXT, format=cast(Any, "xyz"))
    with pytest.raises(ValueError, match="selection"):
        chargefw_io.parse(MOL_TEXT, format="mol", selection="polymers")
    with pytest.raises(ValueError, match="bonds"):
        chargefw_io.parse(MOL_TEXT, format="mol", bonds="explicit")
    with pytest.raises(ValueError):
        chargefw_io.parse(BOND_STRATEGY_PDB, format="pdb", bonds=cast(Any, ""))
    with pytest.raises(ValueError, match="conformers"):
        chargefw_io.parse(MOL_TEXT, format="mol", conformers="first")
