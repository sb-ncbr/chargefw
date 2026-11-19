"""The fields that every issue type has."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TypedDict

import numpy as np
from numpy.typing import NDArray

from ._arrays import index_array
from ._payloads import IssueLocationPayload


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class Issue:
    """A kind, the reason, and where in a molecule the issue is; each issue type narrows ``kind``.

    The message names neither the molecule nor the conformer, nor the method or parameter set.
    """

    kind: str
    message: str
    conformer_index: int | None = None
    # The atom and bond the message names, and the other atom of a pair it names.
    atom_index: int | None = None
    other_atom_index: int | None = None
    bond_index: int | None = None
    # Every atom and bond the issue concerns, in index order.
    atom_indices: NDArray[np.int64] = field(default_factory=index_array)
    bond_indices: NDArray[np.int64] = field(default_factory=index_array)


class LocationArguments(TypedDict):
    conformer_index: int | None
    atom_index: int | None
    other_atom_index: int | None
    bond_index: int | None
    atom_indices: NDArray[np.int64]
    bond_indices: NDArray[np.int64]


def location_arguments(payload: IssueLocationPayload) -> LocationArguments:
    """The location fields of an issue, from its payload."""

    return LocationArguments(
        conformer_index=payload["conformer_index"],
        atom_index=payload["atom_index"],
        other_atom_index=payload["other_atom_index"],
        bond_index=payload["bond_index"],
        atom_indices=index_array(payload["atom_indices"]),
        bond_indices=index_array(payload["bond_indices"]),
    )
