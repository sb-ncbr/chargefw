"""The fitting profiles of SQE, SQE+q0, and SQE+qp, and the change of variables that keeps them stable."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

from ._profile import Profile

if TYPE_CHECKING:
    from .. import Slot

__all__ = ["PROFILES", "StableHardness"]

_SQRT_PI = math.sqrt(math.pi)


@dataclass(frozen=True, slots=True)
class StableHardness:
    """Fit each atom type's hardness as u = √π·hardness·width − 1, with u ≥ 0.

    With every width positive and every bond kappa at or above a positive bound, u ≥ 0 keeps the SQE
    system positive definite for every molecule and geometry: the smallest eigenvalue is then at
    least the smallest kappa. The hardness is calculated back as (1 + u) / (√π·width).
    """

    name: str = "u = √π·hardness·width − 1"
    bounds: tuple[float, float] = (0.0, math.inf)

    def dependencies(self, layout: Sequence[Slot]) -> Mapping[int, tuple[int, ...]]:
        # A set may repeat a key, and each entry has its own hardness and width. The layout lists
        # entries in order, so the n-th hardness of a key belongs with the n-th width of that key.
        widths: dict[object, list[int]] = {}
        for index, slot in enumerate(layout):
            if slot.name == "width":
                widths.setdefault(slot.key, []).append(index)
        seen: dict[object, int] = {}
        pairs = {}
        for index, slot in enumerate(layout):
            if slot.name == "hardness":
                occurrence = seen.get(slot.key, 0)
                seen[slot.key] = occurrence + 1
                if occurrence < len(widths.get(slot.key, ())):
                    pairs[index] = (widths[slot.key][occurrence],)
        return pairs

    def to_internal(
        self, values: NDArray[np.float64], layout: Sequence[Slot]
    ) -> NDArray[np.float64]:
        hardness, width = self._pairs(layout)
        internal = np.array(values, dtype=np.float64)
        internal[hardness] = _SQRT_PI * values[hardness] * values[width] - 1.0
        return internal

    def to_values(
        self, internal: NDArray[np.float64], layout: Sequence[Slot]
    ) -> NDArray[np.float64]:
        hardness, width = self._pairs(layout)
        values = np.array(internal, dtype=np.float64)
        values[hardness] = (1.0 + internal[hardness]) / (_SQRT_PI * internal[width])
        return values

    def derivative(
        self, internal: NDArray[np.float64], layout: Sequence[Slot]
    ) -> NDArray[np.float64]:
        hardness, width = self._pairs(layout)
        result = np.eye(len(layout))
        result[hardness, hardness] = 1.0 / (_SQRT_PI * internal[width])
        result[hardness, width] = -(1.0 + internal[hardness]) / (_SQRT_PI * internal[width] ** 2)
        return result

    def _pairs(self, layout: Sequence[Slot]) -> tuple[list[int], list[int]]:
        pairs = self.dependencies(layout)
        return list(pairs), [width for (width,) in pairs.values()]


_START = {"electronegativity": 0.0, "hardness": 2.0, "width": 1.0, "kappa": 1.0}
# A width must stay positive for the hardness to follow from u; kappa >= 0.01 is the margin that keeps
# the system away from singular.
_BOUNDS = {"width": (0.05, math.inf), "kappa": (0.01, math.inf)}

# The fixed slots remove the invariances: a constant added to every electronegativity, and in SQE+qp
# also to every q0. An atom with one bond moves charge only through that bond, so its hardness and
# the kappa of its bond act as one sum; the hardness of a type whose atoms all have one bond is
# therefore held at its stability limit, u = 0, and the kappa of its bonds take up the rest. When
# those bonds' entries belong to that type alone and their kappa has no upper bound, any stable values
# can be shifted there by raising those kappa, so this costs no accuracy; otherwise nothing is held.
PROFILES = {
    "sqe": Profile(
        bounds=_BOUNDS,
        start=_START,
        fixed=("electronegativity",),
        transforms=(StableHardness(),),
        terminal={"hardness": "kappa"},
        linear=("electronegativity",),
    ),
    "sqeq0": Profile(
        bounds=_BOUNDS,
        start=_START,
        fixed=("electronegativity",),
        transforms=(StableHardness(),),
        terminal={"hardness": "kappa"},
        linear=("electronegativity",),
    ),
    "sqeqp": Profile(
        bounds=_BOUNDS,
        start={**_START, "q0": 0.0},
        fixed=("electronegativity", "q0"),
        transforms=(StableHardness(),),
        terminal={"hardness": "kappa"},
        linear=("electronegativity", "q0"),
    ),
}
