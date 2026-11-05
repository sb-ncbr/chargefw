"""Fitting parameter values of ChargeFW methods to reference charges."""

from __future__ import annotations

import math
import operator
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from typing import SupportsIndex

import numpy as np
from numpy.typing import ArrayLike, NDArray

from ._chargefw import parameters as _native_parameters
from ._methods import Method
from ._parameters import AtomParameterKey, BondParameterKey, ParameterSet, _parameter_set
from ._types import AtomClassifier, BondClassifier, ParameterScope
from .calculation import Plan, _as_collection, parameter_sets
from .core import Molecule, MoleculeCollection

__all__ = [
    "EvaluationError",
    "Evaluator",
    "Slot",
    "make_parameter_template",
]


class EvaluationError(RuntimeError):
    """A calculation with candidate parameter values failed."""


@dataclass(frozen=True, slots=True)
class Slot:
    """One value of a parameter set: a common parameter, or one parameter of one atom or bond type.

    ``key`` is ``None`` for common parameters.
    """

    scope: ParameterScope
    key: AtomParameterKey | BondParameterKey | None
    name: str


def make_parameter_template(
    molecules: Molecule | MoleculeCollection | Iterable[Molecule],
    method: str | Method,
    *,
    atom_types: AtomClassifier,
    bond_types: BondClassifier = "plain",
    id: str = "template",
) -> ParameterSet:
    """Make a parameter set with an entry for every atom and bond type of the molecules.

    Atom types are derived with the ``atom_types`` classifier and bond types with ``bond_types``,
    as described in the parameter-set format. Every entry has the parameter names the method needs,
    with all values zero. Every atom and bond of the molecules matches an entry. The ID must differ
    from the IDs of the bundled parameter sets.
    """

    if id in parameter_sets:
        raise ValueError(f"parameter-set ID {id!r} is the ID of a bundled parameter set")
    method_id = method.id if isinstance(method, Method) else method
    collection = _as_collection(molecules)
    native = _native_parameters._make_parameter_template_catalog(
        collection._native_molecules, method_id, atom_types, bond_types, id
    )
    (descriptor,) = native._descriptors()
    return _parameter_set(descriptor, native)


class Evaluator:
    """Calculates the charges of an assessed plan for any vector of parameter values.

    The plan's molecules are not read, prepared, or matched to parameters again. Values are in the
    order of ``layout``. The plan must use full execution, so that the charges carry no error from
    cutoff or cover execution.
    """

    __slots__ = ("_layout", "_native", "_threads", "_values")

    def __init__(self, plan: Plan, *, threads: int | None = None) -> None:
        parameter_set = plan.parameter_set
        if parameter_set is None:
            raise ValueError(f"method {plan.method.id!r} has no parameters to fit")
        if plan.execution.mode != "full":
            raise ValueError(
                f"fitting needs full execution, but the plan uses {plan.execution.mode}; "
                "assess with execution='full'"
            )
        self._native = plan._native
        self._threads = threads
        self._layout = tuple(
            _slot(scope, index, name, parameter_set)
            for scope, index, name in self._native._layout()
        )
        self._values = np.array(self._native._values(), dtype=np.float64)
        self._values.setflags(write=False)

    @property
    def layout(self) -> tuple[Slot, ...]:
        """The slots, in the order of every value vector."""

        return self._layout

    @property
    def values(self) -> NDArray[np.float64]:
        """The values of the plan's parameter set."""

        return self._values

    def charges(self, values: ArrayLike) -> NDArray[np.float64]:
        """The charges of every assignment, concatenated in assignment order.

        Raises ``EvaluationError`` when a calculation fails, for example with values the method
        cannot use.
        """

        result = self._native._calculate_values(self._vector(values), self._threads)
        if result["status"] != "success" or result["charges"] is None:
            raise EvaluationError(result["failure_message"] or result["status"])
        return result["charges"]

    def jacobian(
        self,
        values: ArrayLike,
        slots: Sequence[SupportsIndex] | None = None,
        bounds: tuple[ArrayLike, ArrayLike] | None = None,
    ) -> NDArray[np.float64]:
        """The derivatives of the charges with respect to the given slots, by forward differences.

        Returns one row per charge and one column per slot, in the order of ``slots`` (default: every
        slot). Slots are indices into ``layout``, from 0; negative indices are not accepted. Every
        step stays within ``bounds``, given as lower and upper values per slot: it goes up when there
        is room, otherwise down, and is shortened to the bound in an interval narrower than the
        step. Each column costs one calculation. Raises ``EvaluationError`` naming the slot when a
        calculation fails.
        """

        point = self._vector(values)
        columns = self._slot_indices(slots)
        lower, upper = self._bounds(bounds)
        if np.any(point < lower) or np.any(point > upper):
            raise ValueError("parameter values must lie within the bounds")
        base = self.charges(point)
        result = np.empty((base.size, len(columns)))
        for column, slot in enumerate(columns):
            moved = point.copy()
            moved[slot] = _step_within_bounds(point[slot], lower[slot], upper[slot])
            step = moved[slot] - point[slot]
            if step == 0.0:
                raise ValueError(
                    f"slot {slot} ({self._describe(slot)}) has no room for a step within its bounds"
                )
            try:
                charges = self.charges(moved)
            except EvaluationError as error:
                raise EvaluationError(
                    f"calculation failed for slot {slot} ({self._describe(slot)}) "
                    f"at value {moved[slot]!r}: {error}"
                ) from error
            result[:, column] = (charges - base) / step
        return result

    def _slot_indices(self, slots: Sequence[SupportsIndex] | None) -> list[int]:
        if slots is None:
            return list(range(len(self._layout)))
        indices = []
        for slot in slots:
            if isinstance(slot, (bool, np.bool_)):
                raise TypeError(f"slot indices must be integers, got {slot!r}")
            try:
                index = operator.index(slot)
            except TypeError:
                raise TypeError(f"slot indices must be integers, got {slot!r}") from None
            if not 0 <= index < len(self._layout):
                raise IndexError(f"slot index {index} is outside 0 to {len(self._layout) - 1}")
            indices.append(index)
        return indices

    def _bounds(
        self, bounds: tuple[ArrayLike, ArrayLike] | None
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        shape = (len(self._layout),)
        if bounds is None:
            return np.full(shape, -math.inf), np.full(shape, math.inf)
        lower = np.broadcast_to(np.asarray(bounds[0], dtype=np.float64), shape)
        upper = np.broadcast_to(np.asarray(bounds[1], dtype=np.float64), shape)
        # Infinite bounds mean no limit; NaN would silently fail every comparison.
        if np.any(np.isnan(lower)) or np.any(np.isnan(upper)):
            raise ValueError("bounds must not be NaN; use infinity for no limit")
        if np.any(lower > upper):
            raise ValueError("lower bounds must not exceed upper bounds")
        return lower, upper

    def _describe(self, slot: int) -> str:
        described = self._layout[slot]
        key = "" if described.key is None else f" of {described.key}"
        return f"{described.scope} {described.name}{key}"

    def _vector(self, values: ArrayLike) -> NDArray[np.float64]:
        vector = np.ascontiguousarray(values, dtype=np.float64)
        if vector.shape != (len(self._layout),):
            raise ValueError(
                f"expected {len(self._layout)} parameter values, got an array of shape "
                f"{vector.shape}"
            )
        if not np.all(np.isfinite(vector)):
            raise ValueError("parameter values must be finite")
        return vector


def _step_within_bounds(value: float, lower: float, upper: float) -> float:
    """The value of a forward-difference step from ``value`` that stays within the bounds."""

    step = math.sqrt(float(np.finfo(np.float64).eps)) * max(1.0, abs(value))
    if upper - value >= step:
        return value + step
    if value - lower >= step:
        return value - step
    return upper if upper - value >= value - lower else lower


def _slot(scope: ParameterScope, index: int, name: str, parameter_set: ParameterSet) -> Slot:
    if scope == "atom":
        return Slot(scope, parameter_set.atom_keys[index], name)
    if scope == "bond":
        return Slot(scope, parameter_set.bond_keys[index], name)
    return Slot(scope, None, name)


Slot.__module__ = "chargefw.parameterisation"
