"""Calculating the charges of an assessed plan for any vector of parameter values."""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass
from typing import SupportsIndex

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .._arrays import as_integer, immutable_array
from .._chargefw import parameterization as _native_parameterization
from .._parameters import AtomParameterKey, BondParameterKey, ParameterSet
from .._types import ParameterScope
from ..calculation import Plan, _native_catalog_for


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


class Evaluator:
    """Calculates the charges of an assessed plan for any vector of parameter values.

    The plan's molecules are not read, prepared, or matched to parameters again. Values are in the
    order of ``layout``. The plan must use full execution, so that the charges carry no error from
    cutoff or cover execution.
    """

    __slots__ = (
        "_assignment_molecules",
        "_assignment_sizes",
        "_atom_entries",
        "_atomic_numbers",
        "_exact",
        "_layout",
        "_matches",
        "_entry_of",
        "_terminal",
        "_native",
        "_parameter_set",
        "_threads",
        "_values",
    )

    def __init__(self, plan: Plan, *, threads: int | None = None) -> None:
        parameter_set = plan.parameter_set
        if parameter_set is None:
            raise ValueError(f"method {plan.method.id!r} has no parameters to fit")
        if plan.execution.mode != "full":
            raise ValueError(
                f"fitting needs full execution, but the plan uses {plan.execution.mode}; "
                "assess with execution='full'"
            )
        self._native = _native_parameterization._NativeEvaluator(plan._native)
        self._parameter_set = parameter_set
        self._threads = threads
        self._layout = tuple(
            _slot(scope, index, name, parameter_set)
            for scope, index, name in self._native._layout()
        )
        self._values = immutable_array(np.asarray(self._native._values(), dtype=np.float64))
        # The method provides exact derivatives, which fit() uses through _normal_equations().
        self._exact = self._native._exact_derivatives()
        # The number of charges and the molecule of each assignment, known after the first
        # successful calculation; they depend only on the molecules and the method, not on values.
        self._assignment_sizes: tuple[int, ...] | None = None
        self._assignment_molecules: tuple[int, ...] | None = None
        self._atomic_numbers = tuple(
            np.array(molecule.atomic_numbers) for molecule in plan._molecules
        )
        self._atom_entries = tuple(np.asarray(entries) for entries in self._native._atom_entries())
        atom_counts, bond_counts = self._native._entry_counts()
        self._matches = tuple(
            atom_counts[index]
            if scope == "atom"
            else bond_counts[index]
            if scope == "bond"
            else None
            for scope, index, _ in self._native._layout()
        )
        # Per atom slot of a type whose atoms all have one bond and whose bonds' entries follow it
        # alone, those bond entries; None otherwise.
        single_bond_entries = self._native._single_bond_entries()
        self._terminal = tuple(
            tuple(single_bond_entries[index][2])
            if scope == "atom"
            and single_bond_entries[index][0]
            and single_bond_entries[index][1]
            and single_bond_entries[index][2]
            else None
            for scope, index, _ in self._native._layout()
        )
        self._entry_of = tuple(index for _, index, _ in self._native._layout())

    @property
    def layout(self) -> tuple[Slot, ...]:
        """The slots, in the order of every value vector."""

        return self._layout

    @property
    def values(self) -> NDArray[np.float64]:
        """The values of the plan's parameter set."""

        return self._values

    @property
    def matches(self) -> tuple[int | None, ...]:
        """Per slot, the number of calculated atoms or bonds of its type; ``None`` for common
        parameters. A type matched by few atoms or bonds is poorly determined by the data."""

        return self._matches

    def charges(self, values: ArrayLike) -> NDArray[np.float64]:
        """The charges of every assignment, concatenated in assignment order.

        Raises ``EvaluationError`` when a calculation fails, for example with values the method
        cannot use.
        """

        result = self._native._calculate_values(self._vector(values), self._threads)
        if (
            result["status"] != "success"
            or result["charges"] is None
            or result["assignment_sizes"] is None
            or result["assignment_molecules"] is None
        ):
            raise EvaluationError(result["failure_message"] or result["status"])
        self._assignment_sizes = tuple(result["assignment_sizes"])
        self._assignment_molecules = tuple(result["assignment_molecules"])
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
        step. If the calculation fails there, the same step the other way is tried, within the
        bounds. Each column costs one calculation, or two after such a failure. Raises
        ``EvaluationError`` naming the slot when no step can be calculated.

        The step is ``sqrt(eps) * max(1, |value|)``, shortened when needed to stay within the bounds,
        so the derivatives are approximate; their accuracy depends on the calculation and the step.
        ``fit()`` uses them only for methods without exact derivatives; for EEM, SQE, SQE+q0, and
        SQE+qp it uses the exact ones.
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
                # Near values a method cannot use, the step the other way may still work.
                first = moved[slot]
                mirrored = point[slot] - step
                if not lower[slot] <= mirrored <= upper[slot] or mirrored == point[slot]:
                    raise EvaluationError(
                        f"calculation failed for slot {slot} ({self._describe(slot)}) "
                        f"at value {moved[slot]!r}: {error}"
                    ) from error
                moved[slot] = mirrored
                step = moved[slot] - point[slot]
                try:
                    charges = self.charges(moved)
                except EvaluationError as mirrored_error:
                    raise EvaluationError(
                        f"calculation failed for slot {slot} ({self._describe(slot)}) at values "
                        f"{first!r} and {mirrored!r}: "
                        f"{mirrored_error}"
                    ) from mirrored_error
            result[:, column] = (charges - base) / step
        return result

    def _normal_equations(
        self,
        values: NDArray[np.float64],
        references: NDArray[np.float64],
        weights: NDArray[np.float64],
    ) -> tuple[NDArray[np.float64], float, NDArray[np.float64], NDArray[np.float64]]:
        """The charges, the weighted sum of squares, JᵀW²J, and JᵀW²r over every slot, from exact
        derivatives; for a method that has them. Raises ``EvaluationError`` when a calculation
        fails."""

        payload = self._native._normal_equations(
            self._vector(values),
            np.ascontiguousarray(references, dtype=np.float64),
            np.ascontiguousarray(weights, dtype=np.float64),
            self._threads,
        )
        if (
            payload["status"] != "success"
            or payload["charges"] is None
            or payload["sum_of_squares"] is None
            or payload["gauss_newton"] is None
            or payload["gradient"] is None
        ):
            raise EvaluationError(payload["failure_message"] or payload["status"])
        size = len(self._layout)
        return (
            payload["charges"],
            float(payload["sum_of_squares"]),
            np.asarray(payload["gauss_newton"]).reshape(size, size),
            np.asarray(payload["gradient"]),
        )

    def _charge_elements(self) -> NDArray[np.int64]:
        """The atomic number of every charge, in the order of ``charges()``."""

        if self._assignment_molecules is None:
            raise RuntimeError("the elements of the charges are known after a calculation")
        return np.concatenate(
            [self._atomic_numbers[molecule] for molecule in self._assignment_molecules]
        )

    def _charge_atom_entries(self) -> NDArray[np.int64]:
        """The atom entry of every charge, in the order of ``charges()``; -1 for none."""

        if self._assignment_molecules is None:
            raise RuntimeError("the atom entries of the charges are known after a calculation")
        return np.concatenate(
            [self._atom_entries[molecule] for molecule in self._assignment_molecules]
        )

    def _slot_indices(self, slots: Sequence[SupportsIndex] | None) -> list[int]:
        if slots is None:
            return list(range(len(self._layout)))
        indices = []
        for slot in slots:
            index = as_integer(slot, "slot index")
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


def _slot_values(
    parameter_set: ParameterSet,
) -> dict[tuple[AtomParameterKey | BondParameterKey | None, str], float]:
    """Every value of a parameter set, by (key, parameter name); the key is None for common
    parameters."""

    catalog = _native_catalog_for(parameter_set)
    layout, values = _native_parameterization._parameter_values(catalog, parameter_set.id)
    return {
        (_slot(scope, index, name, parameter_set).key, name): float(value)
        for (scope, index, name), value in zip(layout, values, strict=True)
    }


def _slot(scope: ParameterScope, index: int, name: str, parameter_set: ParameterSet) -> Slot:
    if scope == "atom":
        return Slot(scope, parameter_set.atom_keys[index], name)
    if scope == "bond":
        return Slot(scope, parameter_set.bond_keys[index], name)
    return Slot(scope, None, name)
