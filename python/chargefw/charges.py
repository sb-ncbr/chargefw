"""Owned charge assignment values."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, SupportsIndex

import numpy as np
from numpy.typing import NDArray

from ._arrays import as_ids, as_integer, immutable_array
from ._types import PortableId
from .core import SourceIdentity


def _charge_values(values: Any) -> NDArray[np.float64]:
    try:
        array = np.asarray(values)
    except Exception as error:
        raise TypeError("values must be a real numeric vector") from error
    if array.ndim != 1:
        raise ValueError("values must be a one-dimensional vector")
    if array.dtype.kind not in "iuf" or array.dtype.kind == "b":
        raise TypeError("values must contain real numbers")
    try:
        normalized = np.asarray(array, dtype=np.float64, order="C")
    except (OverflowError, TypeError, ValueError) as error:
        raise ValueError("values contain values outside the supported range") from error
    if not np.all(np.isfinite(normalized)):
        raise ValueError("values must contain only finite charges")
    return immutable_array(normalized)


def _nonnegative_index(value: SupportsIndex, name: str) -> int:
    result = as_integer(value, name)
    if result < 0:
        raise ValueError(f"{name} must be non-negative")
    return result


@dataclass(frozen=True, slots=True, eq=False)
class ChargeAssignment:
    """Source-mapped, owned charges for one molecule or conformer."""

    values: NDArray[np.float64]
    molecule_index: int
    conformer_index: int | None
    source: SourceIdentity
    atom_ids: tuple[PortableId, ...]

    def __post_init__(self) -> None:
        values = _charge_values(self.values)
        molecule_index = _nonnegative_index(self.molecule_index, "molecule_index")
        conformer_index = (
            None
            if self.conformer_index is None
            else _nonnegative_index(self.conformer_index, "conformer_index")
        )
        if not isinstance(self.source, SourceIdentity):
            raise TypeError("source must be a SourceIdentity")
        atom_ids = as_ids(self.atom_ids, len(values), "atom_ids")

        object.__setattr__(self, "values", values)
        object.__setattr__(self, "molecule_index", molecule_index)
        object.__setattr__(self, "conformer_index", conformer_index)
        object.__setattr__(self, "atom_ids", atom_ids)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ChargeAssignment):
            return NotImplemented
        return (
            np.array_equal(self.values, other.values)
            and self.molecule_index == other.molecule_index
            and self.conformer_index == other.conformer_index
            and self.source == other.source
            and self.atom_ids == other.atom_ids
        )

    def __array__(self, dtype: np.dtype[Any] | None = None, copy: bool | None = None) -> np.ndarray:
        """Expose charge values through NumPy's array conversion protocol."""

        return np.asarray(self.values, dtype=dtype, copy=copy)

    @classmethod
    def _from_native_values(
        cls,
        values: NDArray[np.float64],
        *,
        molecule_index: int,
        conformer_index: int | None,
        source: SourceIdentity,
        atom_ids: tuple[PortableId, ...],
    ) -> ChargeAssignment:
        """Build an assignment from an extension-owned, one-dimensional NumPy array."""
        if values.dtype != np.float64 or values.ndim != 1 or not values.flags.c_contiguous:
            raise RuntimeError("native charge values must be a C-contiguous float64 vector")
        if len(values) != len(atom_ids):
            raise RuntimeError("native charge count does not match the source atom mapping")
        if not np.all(np.isfinite(values)):
            raise RuntimeError("native charge values must be finite")
        if values.flags.writeable:
            raise RuntimeError("native charge values must use immutable storage")
        result = object.__new__(cls)
        object.__setattr__(result, "values", values)
        object.__setattr__(result, "molecule_index", molecule_index)
        object.__setattr__(result, "conformer_index", conformer_index)
        object.__setattr__(result, "source", source)
        object.__setattr__(result, "atom_ids", atom_ids)
        return result
