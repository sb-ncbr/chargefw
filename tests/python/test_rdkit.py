"""Optional RDKit conversion, attachment, and serialization checks."""

from pathlib import Path
from typing import Any, cast, get_type_hints

import chargefw
import numpy as np
import pytest
from chargefw.io import rdkit as chargefw_rdkit

MOL_TEXT = """water
  ChargeFW          3D

  3  2  0  0  0  0  0  0  0  0  1 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    0.9600    0.0000    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.2400    0.9300    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  1  3  1  0  0  0  0
M  END
"""


class FakePoint:
    def __init__(self, x: float, y: float, z: float) -> None:
        self.x = x
        self.y = y
        self.z = z


class FakeConformer:
    def GetId(self) -> int:
        return 0

    def Is3D(self) -> bool:
        return True

    def GetAtomPosition(self, index: int) -> FakePoint:
        return (FakePoint(0.0, 0.0, 0.0), FakePoint(0.96, 0.0, 0.0))[index]


class FakeAtom:
    def __init__(
        self, index: int, atomic_number: int, symbol: str, *, formal_charge: int = 0
    ) -> None:
        self.index = index
        self.atomic_number = atomic_number
        self.symbol = symbol
        self.formal_charge = formal_charge
        self.hydrogen_count = 0
        self.properties: dict[str, float] = {}

    def GetIdx(self) -> int:
        return self.index

    def GetAtomicNum(self) -> int:
        return self.atomic_number

    def GetFormalCharge(self) -> int:
        return self.formal_charge

    def GetSymbol(self) -> str:
        return self.symbol

    def GetTotalNumHs(self) -> int:
        return self.hydrogen_count

    def HasProp(self, name: str) -> bool:
        return name in self.properties

    def GetProp(self, name: str) -> str:
        raise KeyError(name)

    def SetDoubleProp(self, name: str, value: float) -> None:
        self.properties[name] = value


class FakeBondType:
    SINGLE = "SINGLE"
    DOUBLE = "DOUBLE"
    TRIPLE = "TRIPLE"
    AROMATIC = "AROMATIC"
    DATIVEONE = "DATIVEONE"
    DATIVE = "DATIVE"
    DATIVEL = "DATIVEL"
    DATIVER = "DATIVER"
    ZERO = "ZERO"


class FakeBond:
    def __init__(
        self, bond_type: str = FakeBondType.SINGLE, *, aromatic: bool = False, query: bool = False
    ) -> None:
        self.bond_type = bond_type
        self.aromatic = aromatic
        self.query = query

    def GetBeginAtomIdx(self) -> int:
        return 0

    def GetEndAtomIdx(self) -> int:
        return 1

    def HasQuery(self) -> bool:
        return self.query

    def GetBondType(self) -> str:
        return self.bond_type

    def GetIsAromatic(self) -> bool:
        return self.aromatic


class FakeMol:
    def __init__(self, source: "FakeMol | None" = None, quickCopy: bool = False) -> None:
        self.atoms: tuple[FakeAtom, ...] = (
            (FakeAtom(0, 8, "O"), FakeAtom(1, 1, "H")) if source is None else source.atoms
        )
        self.bonds = (FakeBond(),)
        self.conformers: tuple[FakeConformer, ...] = (FakeConformer(),)
        self.properties: dict[str, str] = {}

    def UpdatePropertyCache(self, strict: bool = True) -> None:
        pass

    def GetAtoms(self) -> tuple[FakeAtom, ...]:
        return self.atoms

    def GetBonds(self) -> tuple[FakeBond, ...]:
        return self.bonds

    def GetConformers(self) -> tuple[FakeConformer, ...]:
        return self.conformers

    def HasProp(self, name: str) -> bool:
        return name in self.properties

    def GetProp(self, name: str) -> str:
        return self.properties[name]

    def GetNumAtoms(self) -> int:
        return len(self.atoms)

    def GetAtomWithIdx(self, index: int) -> FakeAtom:
        return self.atoms[index]


class FakeChemistry:
    Mol = FakeMol
    BondType = FakeBondType

    @staticmethod
    def CreateAtomDoublePropertyList(molecule: FakeMol, name: str) -> None:
        molecule.properties[f"atom.dprop.{name}"] = "created"


@pytest.fixture
def fake_rdkit(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(chargefw_rdkit, "_require_rdkit", lambda: FakeChemistry)


@pytest.fixture
def chem() -> Any:
    return pytest.importorskip("rdkit.Chem")


@pytest.mark.usefixtures("fake_rdkit")
def test_conversion_without_conformers_supports_coordinate_independent_methods() -> None:
    target: Any = FakeMol()
    target.conformers = ()

    molecule = chargefw_rdkit.from_mol(target)

    assert molecule.coordinates.shape == (0, 2, 3)
    assert not molecule.has_coordinates
    result = chargefw.calculate(molecule, method="formal")
    np.testing.assert_array_equal(result.assignments[0].values, [0.0, 0.0])


@pytest.mark.usefixtures("fake_rdkit")
def test_conversion_and_charge_attachment() -> None:
    target: Any = FakeMol()
    molecule = chargefw_rdkit.from_mol(target, source_name="water")
    result = chargefw.calculate(molecule, method="formal")
    chargefw_rdkit.attach_charges(target, result)
    with pytest.raises(ValueError):
        chargefw_rdkit.attach_charges(target, result)
    chargefw_rdkit.attach_charges(target, result, overwrite=True)

    assert molecule.atom_ids == (0, 1)
    assert target.atoms[0].properties["ChargeFWPartialCharge"] == 0.0
    assert target.properties["atom.dprop.ChargeFWPartialCharge"] == "created"


@pytest.mark.usefixtures("fake_rdkit")
@pytest.mark.parametrize(
    "bond_type", [FakeBondType.AROMATIC, FakeBondType.DATIVE, FakeBondType.DATIVEONE]
)
def test_conversion_rejects_unsupported_bond_types_by_default(bond_type: str) -> None:
    target: Any = FakeMol()
    target.bonds = (FakeBond(bond_type),)
    with pytest.raises(ValueError, match=bond_type):
        chargefw_rdkit.from_mol(target)


@pytest.mark.usefixtures("fake_rdkit")
@pytest.mark.parametrize(
    "bond_type",
    [
        FakeBondType.AROMATIC,
        FakeBondType.DATIVEONE,
        FakeBondType.DATIVE,
        FakeBondType.DATIVEL,
        FakeBondType.DATIVER,
    ],
)
def test_single_conversion_normalizes_aromatic_and_dative_bonds(bond_type: str) -> None:
    target: Any = FakeMol()
    target.bonds = (FakeBond(bond_type),)
    molecule = chargefw_rdkit.from_mol(target, bond_conversion="single")
    assert molecule.bonds.tolist() == [[0, 1, 1]]


@pytest.mark.usefixtures("fake_rdkit")
def test_single_conversion_rejects_zero_and_query_bonds() -> None:
    target: Any = FakeMol()
    target.bonds = (FakeBond(FakeBondType.ZERO),)
    with pytest.raises(ValueError, match="ZERO"):
        chargefw_rdkit.from_mol(target, bond_conversion="single")

    target.bonds = (FakeBond(query=True),)
    with pytest.raises(ValueError, match="query"):
        chargefw_rdkit.from_mol(target, bond_conversion="single")


@pytest.mark.usefixtures("fake_rdkit")
def test_bond_conversion_is_validated() -> None:
    with pytest.raises(TypeError):
        chargefw_rdkit.from_mol(FakeMol(), bond_conversion=True)  # type: ignore[arg-type]
    with pytest.raises(ValueError):
        chargefw_rdkit.from_mol(FakeMol(), bond_conversion="all")  # type: ignore[arg-type]


@pytest.mark.usefixtures("fake_rdkit")
@pytest.mark.parametrize(
    ("atomic_numbers", "atom_ids"),
    [
        pytest.param([8, 8], [0, 0], id="non-bijective"),
        pytest.param([8, 1], ["0", "1"], id="string"),
    ],
)
def test_attachment_requires_bijective_integer_atom_ids(
    atomic_numbers: list[int], atom_ids: list[Any]
) -> None:
    target: Any = FakeMol()
    target.atoms = tuple(
        FakeAtom(index, number, "O" if number == 8 else "H")
        for index, number in enumerate(atomic_numbers)
    )
    molecule = chargefw.Molecule(atomic_numbers, atom_ids=atom_ids)
    result = chargefw.calculate(molecule, method="formal")

    with pytest.raises(ValueError):
        chargefw_rdkit.attach_charges(target, result)

    assert not any(atom.HasProp("ChargeFWPartialCharge") for atom in target.atoms)


@pytest.mark.usefixtures("fake_rdkit")
def test_attachment_accepts_index_compatible_atom_ids() -> None:
    target: Any = FakeMol()
    target.atoms = (
        FakeAtom(0, 8, "O", formal_charge=-1),
        FakeAtom(1, 1, "H", formal_charge=1),
    )
    molecule = chargefw.Molecule(
        [1, 8],
        formal_charges=[1, -1],
        atom_ids=cast(Any, [np.int64(1), np.int64(0)]),
    )
    result = chargefw.calculate(molecule, method="formal")
    np.testing.assert_array_equal(result.assignments[0].values, [1.0, -1.0])

    chargefw_rdkit.attach_charges(target, result)

    assert target.atoms[0].properties["ChargeFWPartialCharge"] == -1.0
    assert target.atoms[1].properties["ChargeFWPartialCharge"] == 1.0


def test_missing_dependency_is_actionable(monkeypatch: pytest.MonkeyPatch) -> None:
    error = ModuleNotFoundError("No module named 'rdkit'", name="rdkit")

    def missing(name: str) -> Any:
        raise error

    monkeypatch.setattr(chargefw_rdkit, "import_module", missing)
    with pytest.raises(ImportError, match=r"pip install chargefw\[rdkit\]"):
        chargefw_rdkit.from_mol(cast(Any, object()))


def test_public_annotations_resolve_to_runtime_types(chem: Any) -> None:
    hints = get_type_hints(chargefw_rdkit.from_mol)
    assert hints["molecule"] is chem.Mol
    assert hints["return"] is chargefw.Molecule
    hints = get_type_hints(chargefw_rdkit.attach_charges)
    assert hints["molecule"] is chem.Mol
    assert hints["result"] is chargefw.CalculationResult
    assert hints["return"] is type(None)


def test_real_rdkit_conversion_attachment_and_sd_serialization(chem: Any, tmp_path: Path) -> None:
    target = chem.MolFromMolBlock(MOL_TEXT, sanitize=False, removeHs=False)
    assert target is not None
    molecule = chargefw_rdkit.from_mol(target, source_name="water.mol")
    result = chargefw.calculate(molecule, method="formal")

    chargefw_rdkit.attach_charges(target, result)

    assert target.GetAtomWithIdx(0).GetDoubleProp("ChargeFWPartialCharge") == 0.0
    path = tmp_path / "charged.sdf"
    writer = chem.SDWriter(str(path))
    writer.write(target)
    writer.close()
    loaded = chem.SDMolSupplier(str(path), sanitize=False, removeHs=False)[0]
    assert loaded is not None
    assert loaded.GetAtomWithIdx(0).GetDoubleProp("ChargeFWPartialCharge") == 0.0


def test_real_rdkit_conversion_requires_explicit_hydrogens(chem: Any) -> None:
    for molecule in (
        chem.MolFromSmiles("CCO"),
        chem.MolFromSmiles("[NH4+]"),
        chem.MolFromSmiles("CO", sanitize=False),
    ):
        with pytest.raises(ValueError, match="Chem.AddHs"):
            chargefw_rdkit.from_mol(molecule)
    assert chargefw_rdkit.from_mol(chem.AddHs(chem.MolFromSmiles("CCO"))).atom_count == 9


def test_real_rdkit_conversion_refuses_2d_conformers(chem: Any) -> None:
    all_chem = pytest.importorskip("rdkit.Chem.AllChem")
    water = chem.AddHs(chem.MolFromSmiles("O"))
    all_chem.Compute2DCoords(water)
    with pytest.raises(ValueError, match="2D conformer"):
        chargefw_rdkit.from_mol(water)

    # A real flat geometry can be marked as 3D.
    flat = chem.Mol(water)
    flat.GetConformer().Set3D(True)
    assert chargefw_rdkit.from_mol(flat).conformer_count == 1

    # Without conformers, methods that need no geometry still apply.
    without_conformers = chem.Mol(water)
    without_conformers.RemoveAllConformers()
    assert chargefw_rdkit.from_mol(without_conformers).conformer_count == 0
    assert all_chem.EmbedMolecule(water, randomSeed=1) == 0
    assert chargefw_rdkit.from_mol(water).conformer_count == 1


def test_real_rdkit_aromatic_and_dative_bond_conversion(chem: Any) -> None:
    aromatic = chem.AddHs(chem.MolFromSmiles("c1ccccc1"))
    assert aromatic is not None
    with pytest.raises(ValueError, match="AROMATIC"):
        chargefw_rdkit.from_mol(aromatic)
    normalized = chargefw_rdkit.from_mol(aromatic, bond_conversion="single")
    assert np.all(normalized.bonds[:, 2] == 1)

    for name in ("DATIVE", "DATIVEONE"):
        editable = chem.RWMol()
        donor = chem.Atom(7)
        donor.SetNoImplicit(True)
        editable.AddAtom(donor)
        editable.AddAtom(chem.Atom(26))
        editable.AddBond(0, 1, getattr(chem.BondType, name))
        target = editable.GetMol()
        assert target.GetBondWithIdx(0).GetBondTypeAsDouble() == 1.0
        with pytest.raises(ValueError, match=name):
            chargefw_rdkit.from_mol(target)
        normalized = chargefw_rdkit.from_mol(target, bond_conversion="single")
        assert normalized.bonds.tolist() == [[0, 1, 1]]
