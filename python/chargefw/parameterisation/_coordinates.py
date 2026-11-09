"""The profile's changes of variables, applied together."""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence

import numpy as np
from numpy.typing import NDArray

from . import profiles
from ._evaluator import Evaluator, Slot


class _Coordinates:
    """The profile's transforms, applied together: values to internal coordinates and back.

    Transforms act on disjoint slots, so they apply one after another, and the derivative of the
    combination is the product of theirs.
    """

    def __init__(self, transforms: Sequence[profiles.Transform], layout: Sequence[Slot]) -> None:
        self._transforms = tuple(transforms)
        self._layout = layout
        self._owner: dict[int, profiles.Transform] = {}
        self._dependencies: dict[int, tuple[int, ...]] = {}
        for transform in self._transforms:
            for slot, depends_on in transform.dependencies(layout).items():
                if slot in self._owner or any(d in self._owner for d in depends_on):
                    raise ValueError("the transforms of a profile must act on separate slots")
                self._owner[slot] = transform
                self._dependencies[slot] = depends_on

    def check(
        self,
        evaluator: Evaluator,
        free: Sequence[int],
        lower: NDArray[np.float64],
        upper: NDArray[np.float64],
        held: set[int],
    ) -> None:
        """A transformed slot is unbounded, and is free whenever a slot it depends on is free,
        unless the profile holds what the fit varies for it (``held``)."""

        bounded = [s for s in self._owner if math.isfinite(lower[s]) or math.isfinite(upper[s])]
        if bounded:
            raise ValueError(
                "bounds must leave these slots unbounded, because the profile transforms them and "
                "the transform sets their bounds: "
                + "; ".join(f"slot {s} ({self.describe(evaluator)(s)})" for s in bounded)
            )
        free_set = set(free)
        following = [
            s
            for s, depends_on in self._dependencies.items()
            if s not in free_set and s not in held and any(d in free_set for d in depends_on)
        ]
        if following:
            raise ValueError(
                "these slots are not free, but slots their transform depends on are, so their "
                "values would change anyway; free them too, or fix the slots they depend on: "
                + "; ".join(f"slot {s} ({self.describe(evaluator)(s)})" for s in following)
            )

    def bounds(
        self, lower: NDArray[np.float64], upper: NDArray[np.float64]
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        internal_lower, internal_upper = np.array(lower), np.array(upper)
        for slot, transform in self._owner.items():
            internal_lower[slot], internal_upper[slot] = transform.bounds
        return internal_lower, internal_upper

    def to_internal(self, values: NDArray[np.float64]) -> NDArray[np.float64]:
        internal = np.array(values, dtype=np.float64)
        for transform in reversed(self._transforms):
            internal = transform.to_internal(internal, self._layout)
        return internal

    def to_values(self, internal: NDArray[np.float64]) -> NDArray[np.float64]:
        values = np.array(internal, dtype=np.float64)
        for transform in self._transforms:
            values = transform.to_values(values, self._layout)
        return values

    def derivative(self, internal: NDArray[np.float64]) -> NDArray[np.float64]:
        result = np.eye(len(self._layout))
        point = np.array(internal, dtype=np.float64)
        for transform in self._transforms:
            result = transform.derivative(point, self._layout) @ result
            point = transform.to_values(point, self._layout)
        return result

    def plain(self, slot: int) -> bool:
        """Whether the fit varies the slot's value itself, and no transformed slot depends on it."""

        return slot not in self._owner and not any(
            slot in depends_on for depends_on in self._dependencies.values()
        )

    def describe(self, evaluator: Evaluator) -> Callable[[int], str]:
        def described(slot: int) -> str:
            text = evaluator._describe(slot)
            if slot in self._owner:
                text += f", fitted as {self._owner[slot].name}"
            return text

        return described
