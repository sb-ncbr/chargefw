"""Immutable bundled parameter-set descriptors and catalogs."""

from __future__ import annotations

from dataclasses import dataclass

from ._catalog import _Catalog
from ._payloads import (
    AtomParameterKeyPayload,
    BondParameterKeyPayload,
    ParameterSetDescriptorPayload,
)
from ._types import AtomClassifier, BondClassifier


@dataclass(frozen=True, slots=True)
class AtomParameterKey:
    """Atom key of one parameter-set entry."""

    element: str
    classifier: AtomClassifier
    type: str


@dataclass(frozen=True, slots=True)
class BondParameterKey:
    """Bond key of one parameter-set entry; endpoint order is not significant."""

    atoms: tuple[AtomParameterKey, AtomParameterKey]
    classifier: BondClassifier
    type: str


@dataclass(frozen=True, slots=True, repr=False)
class ParameterSet:
    """Metadata for an installed immutable parameter set."""

    id: str
    method: str
    name: str
    publication: str
    notes: str
    atom_keys: tuple[AtomParameterKey, ...]
    bond_keys: tuple[BondParameterKey, ...]

    def __repr__(self) -> str:
        return f"{type(self).__name__}(id={self.id!r}, method={self.method!r}, name={self.name!r})"


class ParameterSetCatalog(_Catalog[ParameterSet]):
    """Installed parameter sets keyed by ID."""

    __slots__ = ()

    def __init__(self, values: tuple[ParameterSet, ...]) -> None:
        super().__init__(values, "parameter-set")


def _atom_key(value: AtomParameterKeyPayload) -> AtomParameterKey:
    return AtomParameterKey(
        element=value["element"], classifier=value["classifier"], type=value["type"]
    )


def _bond_key(value: BondParameterKeyPayload) -> BondParameterKey:
    first, second = value["atoms"]
    return BondParameterKey(
        atoms=(_atom_key(first), _atom_key(second)),
        classifier=value["classifier"],
        type=value["type"],
    )


def _parameter_set(value: ParameterSetDescriptorPayload) -> ParameterSet:
    return ParameterSet(
        id=value["id"],
        method=value["method_id"],
        name=value["name"],
        publication=value["publication"],
        notes=value["notes"],
        atom_keys=tuple(_atom_key(key) for key in value["atom_keys"]),
        bond_keys=tuple(_bond_key(key) for key in value["bond_keys"]),
    )


AtomParameterKey.__module__ = "chargefw"
BondParameterKey.__module__ = "chargefw"
ParameterSet.__module__ = "chargefw"
