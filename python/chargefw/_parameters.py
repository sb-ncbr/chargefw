"""Immutable parameter-set descriptors and catalogs."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from ._catalog import _Catalog
from ._chargefw.parameters import _NativeParameterCatalog
from ._payloads import (
    AtomParameterKeyPayload,
    BondParameterKeyPayload,
    ParameterSetDescriptorPayload,
)
from ._public import set_public_module
from ._types import AtomClassifier, BondClassifier

if TYPE_CHECKING:
    from ._methods import Method


@dataclass(frozen=True, slots=True, kw_only=True)
class AtomParameterKey:
    """Atom key of one parameter-set entry."""

    element: str
    classifier: AtomClassifier
    type: str


@dataclass(frozen=True, slots=True, kw_only=True)
class BondParameterKey:
    """Bond key of one parameter-set entry; endpoint order is not significant."""

    atoms: tuple[AtomParameterKey, AtomParameterKey]
    classifier: BondClassifier
    type: str


class ParameterSet:
    """Metadata for an immutable parameter set, bundled or read from a file."""

    __slots__ = (
        "_id",
        "_method_id",
        "_method",
        "_name",
        "_publication",
        "_notes",
        "_atom_keys",
        "_bond_keys",
        "_native",
    )

    _id: str
    # The ID compares and hashes the set; _link() sets the Method once the methods exist, because
    # each Method lists its bundled sets.
    _method_id: str
    _name: str
    _publication: str
    _notes: str
    _atom_keys: tuple[AtomParameterKey, ...]
    _bond_keys: tuple[BondParameterKey, ...]
    # The native set of a parameter set that is not bundled; None for bundled sets.
    _native: _NativeParameterCatalog | None

    def __init__(self, *args: object, **kwargs: object) -> None:
        raise TypeError(
            "ParameterSet objects come from chargefw.parameter_sets and "
            "chargefw.read_parameter_set()"
        )

    @classmethod
    def _create(
        cls,
        value: ParameterSetDescriptorPayload,
        native: _NativeParameterCatalog | None = None,
    ) -> ParameterSet:
        parameter_set = object.__new__(cls)
        parameter_set._id = value["id"]
        parameter_set._method_id = value["method_id"]
        parameter_set._name = value["name"]
        parameter_set._publication = value["publication"]
        parameter_set._notes = value["notes"]
        parameter_set._atom_keys = tuple(_atom_key(key) for key in value["atom_keys"])
        parameter_set._bond_keys = tuple(_bond_key(key) for key in value["bond_keys"])
        parameter_set._native = native
        return parameter_set

    @property
    def id(self) -> str:
        return self._id

    def _link(self, method: Method) -> None:
        if method.id != self._method_id:
            raise RuntimeError(f"parameter set {self._id!r} is not for method {method.id!r}")
        self._method = method

    @property
    def method(self) -> Method:
        return self._method

    @property
    def name(self) -> str:
        return self._name

    @property
    def publication(self) -> str:
        return self._publication

    @property
    def notes(self) -> str:
        return self._notes

    @property
    def atom_keys(self) -> tuple[AtomParameterKey, ...]:
        return self._atom_keys

    @property
    def bond_keys(self) -> tuple[BondParameterKey, ...]:
        return self._bond_keys

    def _identity(self) -> tuple[object, ...]:
        # The metadata does not include the values, so a set read from a file equals only sets
        # that share its native set, which compares by identity.
        return (
            self._id,
            self._method_id,
            self._name,
            self._publication,
            self._notes,
            self._atom_keys,
            self._bond_keys,
            self._native,
        )

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ParameterSet):
            return NotImplemented
        return self._identity() == other._identity()

    def __hash__(self) -> int:
        return hash(self._identity())

    def __repr__(self) -> str:
        return (
            f"{type(self).__name__}(id={self.id!r}, method={self._method_id!r}, name={self.name!r})"
        )


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


set_public_module("chargefw", AtomParameterKey, BondParameterKey, ParameterSet)
