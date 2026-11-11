"""Immutable bundled periodic-table identities."""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
from operator import index as as_index
from typing import SupportsIndex

from ._chargefw import core as _native_core
from ._payloads import ElementPayload
from ._public import set_public_module


@dataclass(frozen=True, slots=True)
class Element:
    """Identity of one element in the bundled periodic table."""

    atomic_number: int
    symbol: str
    name: str
    period: int
    group: int | None


class PeriodicTable:
    """Bundled elements in atomic-number order, indexed by atomic number or symbol."""

    __slots__ = ("_by_number", "_by_symbol", "_elements")

    _elements: tuple[Element, ...]
    _by_number: dict[int, Element]
    _by_symbol: dict[str, Element]

    def __init__(self, elements: tuple[Element, ...]) -> None:
        object.__setattr__(self, "_elements", elements)
        object.__setattr__(
            self, "_by_number", {element.atomic_number: element for element in elements}
        )
        object.__setattr__(self, "_by_symbol", {element.symbol: element for element in elements})

    def __setattr__(self, name: str, value: object) -> None:
        raise AttributeError(f"{type(self).__name__} is immutable")

    def __len__(self) -> int:
        return len(self._elements)

    def __iter__(self) -> Iterator[Element]:
        return iter(self._elements)

    def __getitem__(self, key: SupportsIndex | str) -> Element:
        """Return an element by atomic number, including NumPy integers, or by exact symbol."""

        element = self._lookup(key)
        if element is None:
            raise KeyError(f"unknown element: {key!r}")
        return element

    def __contains__(self, key: object) -> bool:
        try:
            return self._lookup(key) is not None
        except TypeError:
            return False

    def _lookup(self, key: object) -> Element | None:
        if isinstance(key, str):
            return self._by_symbol.get(key)
        if isinstance(key, bool) or not isinstance(key, SupportsIndex):
            raise TypeError("element keys must be atomic numbers or symbols")
        return self._by_number.get(as_index(key))

    def __repr__(self) -> str:
        return f"{type(self).__name__}(size={len(self)})"


def _element(value: ElementPayload) -> Element:
    return Element(
        atomic_number=value["atomic_number"],
        symbol=value["symbol"],
        name=value["name"],
        period=value["period"],
        group=value["group"],
    )


periodic_table = PeriodicTable(tuple(_element(value) for value in _native_core._elements()))

set_public_module("chargefw", Element, PeriodicTable)
