"""Fitting parameter values of ChargeFW methods to reference charges."""

from __future__ import annotations

import math
import operator
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import SupportsIndex

import numpy as np
from numpy.typing import ArrayLike, NDArray

from ._chargefw import parameters as _native_parameters
from ._methods import Method
from ._parameters import AtomParameterKey, BondParameterKey, ParameterSet, _parameter_set
from ._types import AtomClassifier, BondClassifier, ParameterScope
from .calculation import Plan, _as_collection, _default_parameter_catalog, parameter_sets
from .core import Molecule, MoleculeCollection

__all__ = [
    "EvaluationError",
    "Evaluator",
    "FitResult",
    "Slot",
    "fit",
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

    __slots__ = (
        "_assignment_sizes",
        "_layout",
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
        self._native = plan._native
        self._parameter_set = parameter_set
        self._threads = threads
        self._layout = tuple(
            _slot(scope, index, name, parameter_set)
            for scope, index, name in self._native._layout()
        )
        self._values = np.array(self._native._values(), dtype=np.float64)
        self._values.setflags(write=False)
        # The number of charges in each assignment, known after the first successful calculation;
        # it depends only on the molecules and the method, not on the values.
        self._assignment_sizes: tuple[int, ...] | None = None

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
        if (
            result["status"] != "success"
            or result["charges"] is None
            or result["assignment_sizes"] is None
        ):
            raise EvaluationError(result["failure_message"] or result["status"])
        self._assignment_sizes = tuple(result["assignment_sizes"])
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


@dataclass(frozen=True, slots=True)
class FitResult:
    """The outcome of ``fit()``.

    ``values`` and ``start`` hold every slot in the order of the evaluator's layout; ``free`` lists
    the slots the fit changed. ``rmsd`` compares the charges at ``values`` with the references,
    without weights. ``failed_evaluations`` counts trial values whose calculation failed and were
    therefore rejected. ``parameter_set()`` turns the fitted values into a parameter set.
    """

    values: NDArray[np.float64]
    start: NDArray[np.float64]
    free: tuple[int, ...]
    success: bool
    message: str
    evaluations: int
    failed_evaluations: int
    rmsd: float
    # The parameter set the fit started from, with its original values; parameter_set() combines its
    # keys with ``values``. It keeps no training data.
    _parameter_set: ParameterSet = field(repr=False, compare=False)

    def parameter_set(
        self, *, id: str, name: str = "", publication: str = "", notes: str = ""
    ) -> ParameterSet:
        """The fitted values as a parameter set for the evaluator's method.

        Pass it to ``assess()`` or ``calculate()``, or write it with ``write_parameter_set()``.
        The ID must differ from the IDs of the bundled parameter sets.
        """

        if id in parameter_sets:
            raise ValueError(f"parameter-set ID {id!r} is the ID of a bundled parameter set")
        source = self._parameter_set
        catalog = source._native or _default_parameter_catalog()
        native = catalog._with_values(
            source.id, np.ascontiguousarray(self.values), id, name, publication, notes
        )
        (descriptor,) = native._descriptors()
        return _parameter_set(descriptor, native)


def fit(
    evaluator: Evaluator,
    references: Sequence[ArrayLike] | ArrayLike,
    *,
    weights: Sequence[ArrayLike] | ArrayLike | None = None,
    start: ArrayLike
    | Mapping[tuple[AtomParameterKey | BondParameterKey | None, str], float]
    | None = None,
    free: Sequence[SupportsIndex] | None = None,
    bounds: tuple[ArrayLike, ArrayLike] | None = None,
    max_evaluations: int | None = None,
) -> FitResult:
    """Fit parameter values so that the charges match reference charges in the least-squares sense.

    ``references`` holds the reference charges of every assignment, in assignment order: one array
    per assignment, or one array for all. ``weights`` has the same shape, or is omitted for equal
    weights. ``start`` is a vector in layout order, or a mapping from (key, parameter name) to the
    values that differ from the evaluator's values; the key is ``None`` for common parameters.
    ``free`` lists the slots to fit (default: all); the others keep their start values. ``bounds``
    gives lower and upper values per slot.

    One array per assignment must match the assignments in number and size. ``free`` must not
    repeat a slot.

    Before fitting, the start is checked: it must lie strictly inside the bounds of every free slot,
    at least ``1e-10 * max(1, |bound|)`` from a finite bound, so that SciPy uses it unchanged; the
    charges must be calculable; and every free slot must change them. Calculations that fail at
    trial values are rejected, and the optimizer tries a shorter step. A Jacobian step that fails is
    tried the other way; if that fails too, the fit stops with ``EvaluationError``. Needs SciPy.
    """

    try:
        from scipy.optimize import least_squares  # type: ignore[import-untyped, unused-ignore]
    except ImportError as error:
        raise ImportError(
            "fit() needs SciPy; install it, for example with pip install chargefw[parameterisation]"
        ) from error

    free_slots = evaluator._slot_indices(free)
    if not free_slots:
        raise ValueError("fit needs at least one free slot")
    if len(set(free_slots)) != len(free_slots):
        raise ValueError("free slots must not repeat")
    lower, upper = evaluator._bounds(bounds)
    initial = _start_values(evaluator, start)
    _check_start_in_bounds(evaluator, initial, free_slots, lower, upper)

    try:
        start_charges = evaluator.charges(initial)
    except EvaluationError as error:
        raise EvaluationError(f"the charges cannot be calculated at the start: {error}") from error
    assert evaluator._assignment_sizes is not None
    sizes = evaluator._assignment_sizes
    reference = _flat(references, "references", sizes)
    weight = np.ones(start_charges.size) if weights is None else _flat(weights, "weights", sizes)

    start_jacobian = evaluator.jacobian(initial, free_slots, (lower, upper))
    insensitive = [
        slot for column, slot in enumerate(free_slots) if not start_jacobian[:, column].any()
    ]
    if insensitive:
        raise ValueError(
            "the charges do not change with these free slots at the start, so the fit cannot "
            "determine them: "
            + "; ".join(f"slot {slot} ({evaluator._describe(slot)})" for slot in insensitive)
            + ". Give a start where they matter, for example different values for different types, "
            "or do not fit them"
        )

    failed = 0

    def values_of(free_values: NDArray[np.float64]) -> NDArray[np.float64]:
        values = initial.copy()
        values[free_slots] = free_values
        return values

    def residuals(free_values: NDArray[np.float64]) -> NDArray[np.float64]:
        nonlocal failed
        try:
            charges = evaluator.charges(values_of(free_values))
        except EvaluationError:
            # Non-finite residuals make the TRF method reject the trial and shorten its step.
            failed += 1
            return np.full(reference.size, np.nan)
        return weight * (charges - reference)

    def jacobian(free_values: NDArray[np.float64]) -> NDArray[np.float64]:
        values = values_of(free_values)
        if np.array_equal(values, initial):
            charges_jacobian = start_jacobian
        else:
            charges_jacobian = evaluator.jacobian(values, free_slots, (lower, upper))
        return weight[:, np.newaxis] * charges_jacobian

    solution = least_squares(
        residuals,
        initial[free_slots],
        jac=jacobian,
        bounds=(lower[free_slots], upper[free_slots]),
        method="trf",
        x_scale="jac",
        max_nfev=max_evaluations,
    )
    fitted = values_of(solution.x)
    fitted.setflags(write=False)
    initial.setflags(write=False)
    final_charges = evaluator.charges(fitted)
    return FitResult(
        values=fitted,
        start=initial,
        free=tuple(free_slots),
        success=bool(solution.success),
        message=str(solution.message),
        evaluations=int(solution.nfev),
        failed_evaluations=failed,
        rmsd=float(np.sqrt(np.mean((final_charges - reference) ** 2))),
        _parameter_set=evaluator._parameter_set,
    )


def _start_values(
    evaluator: Evaluator,
    start: ArrayLike
    | Mapping[tuple[AtomParameterKey | BondParameterKey | None, str], float]
    | None,
) -> NDArray[np.float64]:
    if start is None:
        return np.array(evaluator.values)
    if isinstance(start, Mapping):
        values = np.array(evaluator.values)
        positions = {(slot.key, slot.name): index for index, slot in enumerate(evaluator.layout)}
        for (key, name), value in start.items():
            if (key, name) not in positions:
                raise KeyError(f"no slot for parameter {name!r} of {key!r} in the layout")
            values[positions[(key, name)]] = value
        return evaluator._vector(values).copy()
    return evaluator._vector(start).copy()


def _check_start_in_bounds(
    evaluator: Evaluator,
    start: NDArray[np.float64],
    free: Sequence[int],
    lower: NDArray[np.float64],
    upper: NDArray[np.float64],
) -> None:
    # SciPy's TRF method moves a start closer to a finite bound than 1e-10 * max(1, |bound|) away
    # from it. A start at least that far inside is used exactly as given.
    def margin(bound: float) -> float:
        return 1e-10 * max(1.0, abs(bound)) if math.isfinite(bound) else 0.0

    outside = [
        slot
        for slot in free
        if not lower[slot] + margin(lower[slot]) < start[slot] < upper[slot] - margin(upper[slot])
    ]
    if outside:
        raise ValueError(
            "the start must lie strictly inside the bounds of every free slot, at least "
            "1e-10 * max(1, |bound|) from a finite bound; change the start or the bounds of: "
            + "; ".join(
                f"slot {slot} ({evaluator._describe(slot)}) = {start[slot]!r} in "
                f"[{lower[slot]!r}, {upper[slot]!r}]"
                for slot in outside
            )
        )


def _flat(
    value: Sequence[ArrayLike] | ArrayLike, name: str, sizes: Sequence[int]
) -> NDArray[np.float64]:
    """One value per calculated charge, from one array per assignment or one array for all."""

    if (
        isinstance(value, Sequence)
        and not isinstance(value, str)
        and any(np.ndim(part) > 0 for part in value)
    ):
        if len(value) != len(sizes):
            raise ValueError(
                f"{name} has {len(value)} arrays, but the calculation has {len(sizes)} assignments"
            )
        parts = [np.asarray(part, dtype=np.float64) for part in value]
        for index, (part, size) in enumerate(zip(parts, sizes, strict=True)):
            if part.shape != (size,):
                raise ValueError(
                    f"{name} for assignment {index} must have {size} values, one per atom, not "
                    f"an array of shape {part.shape}"
                )
        flat = np.concatenate(parts)
    else:
        flat = np.ravel(np.asarray(value, dtype=np.float64))
        if flat.size != sum(sizes):
            raise ValueError(
                f"{name} must have {sum(sizes)} values, one per calculated charge, not {flat.size}"
            )
    if not np.all(np.isfinite(flat)):
        raise ValueError(f"{name} must be finite")
    return flat


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
