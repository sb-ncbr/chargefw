"""Owned molecule and collection model checks."""

from dataclasses import FrozenInstanceError
from gc import collect
from typing import Any, cast

import chargefw
import numpy as np
import pytest
from chargefw._chargefw import core as _native_core


def test_imported_source_mapping_is_owned_and_read_only() -> None:
    contents = """\
@<TRIPOS>MOLECULE
mapping
2 0 0 0 0
SMALL
NO_CHARGES
@<TRIPOS>ATOM
001 C1 0 0 0 C.3
10 O1 1 0 0 O.2
@<TRIPOS>BOND
"""
    collection = chargefw.io.parse(contents, format="mol2", source_name="mapping.mol2")
    molecule = collection[0]
    del collection, contents
    collect()

    assert molecule.atom_ids == ("001", "10")
    mapping = molecule.source_mapping
    assert mapping is not None
    assert mapping.format == "mol2"
    assert mapping.source_connectivity == "explicitly-empty"
    assert mapping.atoms == (
        chargefw.SourceAtomReference(0, "001"),
        chargefw.SourceAtomReference(1, "10"),
    )
    assert mapping.conformers[0].sites == mapping.atoms
    assert mapping.components == ()
    with pytest.raises(FrozenInstanceError):
        setattr(mapping, "format", "mol")

    manual = chargefw.Molecule(
        molecule.atomic_numbers,
        coordinates=molecule.coordinates,
        atom_ids=molecule.atom_ids,
    )
    assert manual.source_mapping is None


def test_private_binding_rejects_native_integer_overflow() -> None:
    with pytest.raises(ValueError):
        _native_core._make_molecule(
            np.array([2**40], dtype=np.int64),
            np.zeros(1, dtype=np.int64),
            np.empty((0, 3), dtype=np.int64),
            np.empty((0, 1, 3), dtype=np.float64),
            [""],
            [],
            "",
        )


def test_molecule_owns_normalized_input() -> None:
    atomic_source = np.array([8, 1, 1, 6], dtype=np.int8)
    formal = np.array([0, 0, 0, -1], dtype=np.int8)
    bonds = np.array([[0, 1, 1], [0, 2, 1], [0, 3, 2]], dtype=np.int32)
    coordinates_source = np.arange(24, dtype=np.float32).reshape(2, 4, 3)
    atom_ids = ["O", 17, np.int64(23), "H"]

    molecule = chargefw.Molecule(
        atomic_source,
        formal_charges=formal,
        bonds=bonds,
        coordinates=coordinates_source[:, :, ::-1],
        name="water-like",
        atom_names=["O", "H1", "H2", "C"],
        conformer_names=["a", "b"],
        source_name="fixture.sdf",
        record_index=3,
        record_id="record-4",
        atom_ids=cast(Any, atom_ids),
    )

    atomic_source[:] = 1
    formal[:] = 9
    bonds[:] = 0
    coordinates_source[:] = 99

    np.testing.assert_array_equal(molecule.atomic_numbers, [8, 1, 1, 6])
    np.testing.assert_array_equal(molecule.formal_charges, [0, 0, 0, -1])
    np.testing.assert_array_equal(molecule.bonds, [[0, 1, 1], [0, 2, 1], [0, 3, 2]])
    assert molecule.coordinates.shape == (2, 4, 3)
    assert molecule.coordinates[0, 0].tolist() == [2.0, 1.0, 0.0]
    assert molecule.source == chargefw.SourceIdentity("fixture.sdf", 3, "record-4")
    assert molecule.atom_ids == ("O", 17, 23, "H")
    assert len(molecule) == 4

    readonly = molecule.atomic_numbers
    # Each access returns a fresh view, so changing its metadata leaves the molecule unchanged.
    readonly.shape = (2, 2)
    assert molecule.atomic_numbers.shape == (4,)
    with pytest.raises(ValueError):
        readonly[0] = 2
    with pytest.raises(ValueError):
        readonly.setflags(write=True)
    with pytest.raises(ValueError):
        molecule.formal_charges[0] = 2
    with pytest.raises(ValueError):
        molecule.bonds[0, 0] = 1
    with pytest.raises(ValueError):
        molecule.coordinates[0, 0, 0] = 1.0
    with pytest.raises(AttributeError):
        setattr(molecule, "name", "changed")


def test_collection_is_an_immutable_sequence() -> None:
    molecule = chargefw.Molecule([1])
    collection = chargefw.MoleculeCollection([molecule], name="fixture")

    assert len(collection) == 1
    assert collection[0] is molecule
    assert collection[:] == (molecule,)
    assert tuple(collection) == (molecule,)
    assert collection.molecules == (molecule,)
    assert collection.name == "fixture"


def test_coordinate_defaults_and_empty_molecule() -> None:
    no_coordinates = chargefw.Molecule([1])
    assert no_coordinates.coordinates.shape == (0, 1, 3)
    assert no_coordinates.conformer_count == 0
    assert not no_coordinates.has_coordinates
    assert chargefw.Molecule([]).atomic_numbers.dtype == np.dtype(np.int64)

    one_conformer = chargefw.Molecule([1], coordinates=[[0.0, 0.0, 0.0]])
    assert one_conformer.coordinates.shape == (1, 1, 3)
    assert one_conformer.conformer_count == 1


def test_integer_iterables_and_empty_bonds_are_normalized() -> None:
    molecule = chargefw.Molecule((atomic_number for atomic_number in [8, 1, 1]), bonds=[])
    np.testing.assert_array_equal(molecule.atomic_numbers, [8, 1, 1])
    assert molecule.bonds.shape == (0, 3)


def test_object_integer_arrays_are_range_checked() -> None:
    molecule = chargefw.Molecule(np.array([1, 8], dtype=object))
    np.testing.assert_array_equal(molecule.atomic_numbers, [1, 8])
    with pytest.raises(ValueError):
        chargefw.Molecule(np.array([2**100], dtype=object))


def test_supported_atomic_number_boundaries() -> None:
    assert chargefw.Molecule([100]).atomic_numbers[0] == 100
    with pytest.raises(ValueError, match=r"range 1\.\.100"):
        chargefw.Molecule([101])


def test_periodic_table_identifies_supported_elements() -> None:
    table = chargefw.periodic_table
    assert [element.atomic_number for element in table] == list(range(1, 101))
    molecule = chargefw.Molecule([6, 34])
    assert [table[number].symbol for number in molecule.atomic_numbers] == ["C", "Se"]
    assert table["C"] is table[6]
    assert table["C"].group == 14
    assert table["Ce"].group is None
    assert 0 not in table
    assert "c" not in table
    with pytest.raises(KeyError):
        table[101]
    with pytest.raises(TypeError):
        table[cast(Any, True)]
    with pytest.raises(AttributeError):
        cast(Any, table).extra = None


def test_invalid_inputs_are_rejected() -> None:
    invalid_cases = (
        (TypeError, lambda: chargefw.Molecule([1.0])),
        (ValueError, lambda: chargefw.Molecule([0])),
        (ValueError, lambda: chargefw.Molecule([1, 1], bonds=[[0, 1, 4]])),
        (
            ValueError,
            lambda: chargefw.Molecule([1, 1], bonds=[[0, 1, 1], [1, 0, 1]]),
        ),
        (
            ValueError,
            lambda: chargefw.Molecule([1], coordinates=np.empty((0, 1, 2))),
        ),
        (TypeError, lambda: chargefw.Molecule([1], atom_ids=cast(Any, [[]]))),
        (TypeError, lambda: chargefw.Molecule([1], atom_ids=cast(Any, [True]))),
        (ValueError, lambda: chargefw.Molecule([1], atom_ids=[2**100])),
        (TypeError, lambda: chargefw.Molecule([1], atom_names="H")),
        (ValueError, lambda: chargefw.SourceIdentity(record_index=-1)),
        (TypeError, lambda: chargefw.SourceIdentity(record_index=cast(Any, 1.5))),
        (TypeError, lambda: chargefw.SourceIdentity(record_index=True)),
        (TypeError, lambda: chargefw.SourceIdentity(record_id=cast(Any, ("record", 1)))),
        (
            TypeError,
            lambda: chargefw.Molecule([1], record_index=cast(Any, np.bool(True))),
        ),
    )
    for error_type, operation in invalid_cases:
        with pytest.raises(error_type):
            operation()
