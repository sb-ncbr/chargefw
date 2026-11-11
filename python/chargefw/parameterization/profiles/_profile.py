"""What a profile is, independent of any method."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import TYPE_CHECKING, Protocol

from numpy.typing import NDArray

if TYPE_CHECKING:
    import numpy as np

    from .. import Slot


class Transform(Protocol):
    """A change of variables: for some slots, the fit varies another quantity than the value itself.

    ``name`` says what the fit varies, for messages. ``bounds`` gives its lower and upper value, the same
    for every transformed slot. ``dependencies`` maps each transformed slot to the other slots its value
    depends on. ``to_internal`` and ``to_values`` convert a vector of every slot's value into the fitted
    variables and back, and ``derivative`` gives d values / d internal at an internal vector, a square
    matrix in layout order.
    """

    @property
    def name(self) -> str: ...

    @property
    def bounds(self) -> tuple[float, float]: ...

    def dependencies(self, layout: Sequence[Slot]) -> Mapping[int, tuple[int, ...]]: ...

    def to_internal(
        self, values: NDArray[np.float64], layout: Sequence[Slot]
    ) -> NDArray[np.float64]: ...

    def to_values(
        self, internal: NDArray[np.float64], layout: Sequence[Slot]
    ) -> NDArray[np.float64]: ...

    def derivative(
        self, internal: NDArray[np.float64], layout: Sequence[Slot]
    ) -> NDArray[np.float64]: ...


@dataclass(frozen=True, slots=True)
class Profile:
    """What a fit uses for one method unless told otherwise.

    ``bounds`` maps a parameter name to its lower and upper value, for every slot of that name.
    ``start`` maps a parameter name to a starting value, used for a parameter whose values are all
    zero in the plan's parameter set, as in a template. ``fixed`` names parameters whose first slot
    is held at its start value: one per invariance of the method, a change of values that changes no
    charge. ``transforms`` are changes of variables the fit applies; a transformed slot takes its
    bounds from the transform, not from ``bounds``. ``terminal`` maps an atom parameter to a bond
    parameter that acts on the charges only in its sum with it at an atom with one bond, as SQE's
    hardness and kappa do. For an atom type whose atoms all have one bond, and whose bonds' entries
    belong to no other type, the atom parameter is held at the lower bound of what the fit varies for
    it, and the bond parameter of those bonds takes up the difference: this removes the change of
    both that changes no charge. ``linear`` names parameters on which the charges depend linearly
    when all other values are fixed, such as electronegativities. Unless told which slots to fit,
    the fit first sets their free, unbounded slots to the values that fit the references best for
    the start of the others, and then fits them with the rest. They need no start; only
    combinations of them that the data does not determine keep the start they had.
    """

    bounds: Mapping[str, tuple[float, float]] = field(default_factory=dict)
    start: Mapping[str, float] = field(default_factory=dict)
    fixed: tuple[str, ...] = ()
    transforms: tuple[Transform, ...] = ()
    terminal: Mapping[str, str] = field(default_factory=dict)
    linear: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(self, "bounds", MappingProxyType(dict(self.bounds)))
        object.__setattr__(self, "start", MappingProxyType(dict(self.start)))
        object.__setattr__(self, "terminal", MappingProxyType(dict(self.terminal)))
