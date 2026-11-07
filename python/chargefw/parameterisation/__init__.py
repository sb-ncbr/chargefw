"""Fitting parameter values of ChargeFW methods to reference charges."""

from __future__ import annotations

import math
import operator
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Literal, SupportsIndex

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .._chargefw import parameters as _native_parameters
from .._elements import periodic_table
from .._methods import Method
from .._parameters import AtomParameterKey, BondParameterKey, ParameterSet, _parameter_set
from .._types import AtomClassifier, BondClassifier, ParameterScope
from ..calculation import Plan, _as_collection, _default_parameter_catalog, parameter_sets
from ..core import Molecule, MoleculeCollection
from . import profiles

__all__ = [
    "Diagnostics",
    "ElementDiagnostics",
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
        "_assignment_molecules",
        "_assignment_sizes",
        "_atomic_numbers",
        "_layout",
        "_matches",
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
        # The number of charges and the molecule of each assignment, known after the first
        # successful calculation; they depend only on the molecules and the method, not on values.
        self._assignment_sizes: tuple[int, ...] | None = None
        self._assignment_molecules: tuple[int, ...] | None = None
        self._atomic_numbers = tuple(
            np.array(molecule.atomic_numbers) for molecule in plan._molecules
        )
        atom_counts, bond_counts = self._native._entry_counts()
        self._matches = tuple(
            atom_counts[index]
            if scope == "atom"
            else bond_counts[index]
            if scope == "bond"
            else None
            for scope, index, _ in self._native._layout()
        )

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

    def _charge_elements(self) -> NDArray[np.int64]:
        """The atomic number of every charge, in the order of ``charges()``."""

        if self._assignment_molecules is None:
            raise RuntimeError("the elements of the charges are known after a calculation")
        return np.concatenate(
            [self._atomic_numbers[molecule] for molecule in self._assignment_molecules]
        )

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
class ElementDiagnostics:
    """How well the charges of one element match the references.

    ``r2`` is ``nan`` when every reference of the element is the same, which leaves it undefined.
    """

    count: int
    rmsd: float
    r2: float


@dataclass(frozen=True, slots=True)
class Diagnostics:
    """What shows whether the fitted values can be trusted.

    ``r2`` is the coefficient of determination of the charges at the fitted values against the
    references, without weights; ``nan`` when every reference is the same. ``at_bounds`` lists free
    slots whose fitted value is at a bound, within ``1e-8 * max(1, |bound|)``; for a slot the profile
    transforms, the bound of the transformed quantity counts, such as the stability limit of the SQE
    hardness. ``failed_messages`` holds the distinct messages of trial calculations that failed, at
    most 20.

    ``singular_values`` come from the weighted Jacobian at the fitted values, for the free slots,
    with each column scaled to unit length, so they compare slots whatever their units; they are in
    decreasing order. ``rank`` counts those above ``1e-6`` times the largest, well above the noise of
    a finite-difference Jacobian. A rank below the number of free slots means the data cannot
    determine some combination of them. ``weakest`` lists the
    free slots that dominate the direction of the smallest singular value: those whose share of it is
    at least half the largest share.

    ``matches`` gives, for each free slot in the order of ``FitResult.free``, the number of
    calculated atoms or bonds of its type, or ``None`` for a common parameter. ``by_element`` maps
    each element symbol, in atomic-number order, to the number of charges of that element and their
    RMSD and R² against the references.
    """

    r2: float
    by_element: Mapping[str, ElementDiagnostics]
    matches: tuple[int | None, ...]
    at_bounds: tuple[int, ...]
    failed_messages: tuple[str, ...]
    singular_values: NDArray[np.float64]
    rank: int
    weakest: tuple[int, ...]


@dataclass(frozen=True, slots=True)
class FitResult:
    """The outcome of ``fit()``.

    ``values`` and ``start`` hold every slot in the order of the evaluator's layout; ``free`` lists
    the slots the fit changed, and ``bounds`` the lower and upper values it used per slot; a slot the
    profile transforms has no bounds of its own there, as its bound applies to the transformed
    quantity. ``rmsd``
    compares the charges at ``values`` with the references, without weights.
    ``failed_evaluations`` counts trial values whose calculation failed and were therefore rejected.
    ``diagnostics`` shows whether the values can be trusted. ``parameter_set()`` turns the fitted
    values into a parameter set.
    """

    values: NDArray[np.float64]
    start: NDArray[np.float64]
    free: tuple[int, ...]
    bounds: tuple[NDArray[np.float64], NDArray[np.float64]]
    success: bool
    message: str
    evaluations: int
    failed_evaluations: int
    rmsd: float
    diagnostics: Diagnostics
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
    profile: profiles.Profile | Literal["default"] | None = "default",
    max_evaluations: int | None = None,
) -> FitResult:
    """Fit parameter values so that the charges match reference charges in the least-squares sense.

    ``references`` holds the reference charges of every assignment, in assignment order: one array
    per assignment, or one array for all. ``weights`` has the same shape, or is omitted for equal
    weights. ``start`` is a vector in layout order, or a mapping from (key, parameter name) to the
    values that differ from the evaluator's values; the key is ``None`` for common parameters.
    ``free`` lists the slots to fit; the others keep their start values. ``bounds`` gives lower and
    upper values per slot.

    ``profile`` gives the defaults for what is omitted: by default the method's entry in
    ``profiles.PROFILES``, or nothing for a method without one, or for ``profile=None``. Without ``free``,
    every slot is fitted except the profile's fixed ones. Without ``bounds``, the profile's bounds
    apply, and no bounds without a profile. A parameter whose values are all zero in the plan's set,
    as in a template, starts at the profile's value for it; a free slot of such a parameter without
    a profile value needs a value in ``start``.

    The profile's transforms change what the fit varies for some slots, such as the SQE hardness,
    fitted as u = √π·hardness·width − 1 ≥ 0. A transformed slot takes its bounds from the transform:
    ``bounds`` must leave it unbounded, and if the slot is not free, the slots it depends on (the
    width) must not be free either.

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

    chosen = (
        profiles.PROFILES.get(evaluator._parameter_set.method) if profile == "default" else profile
    )
    layout = evaluator.layout
    if free is not None:
        free_slots = evaluator._slot_indices(free)
    else:
        fixed = _fixed_slots(layout, chosen)
        free_slots = [slot for slot in range(len(layout)) if slot not in fixed]
    if not free_slots:
        raise ValueError("fit needs at least one free slot")
    if len(set(free_slots)) != len(free_slots):
        raise ValueError("free slots must not repeat")
    lower, upper = evaluator._bounds(
        bounds if bounds is not None else _profile_bounds(layout, chosen)
    )
    coordinates = _Coordinates(chosen.transforms if chosen is not None else (), layout)
    coordinates.check(evaluator, free_slots, lower, upper)
    # The fit works in internal coordinates: the values, except that transformed slots hold the
    # transformed quantity, which carries the transform's bounds.
    internal_lower, internal_upper = coordinates.bounds(lower, upper)
    initial, missing = _start_values(evaluator, start, chosen)
    missing_free = [slot for slot in free_slots if slot in missing]
    if missing_free:
        raise ValueError(
            "these free slots have no start value, because their parameter is not given in the "
            "parameter set and has no default; give them in start: "
            + "; ".join(f"slot {slot} ({evaluator._describe(slot)})" for slot in missing_free)
        )
    internal_initial = coordinates.to_internal(initial)
    _check_start_in_bounds(
        coordinates.describe(evaluator),
        internal_initial,
        free_slots,
        internal_lower,
        internal_upper,
    )

    try:
        start_charges = evaluator.charges(initial)
    except EvaluationError as error:
        raise EvaluationError(f"the charges cannot be calculated at the start: {error}") from error
    assert evaluator._assignment_sizes is not None
    sizes = evaluator._assignment_sizes
    reference = _flat(references, "references", sizes)
    weight = np.ones(start_charges.size) if weights is None else _flat(weights, "weights", sizes)

    start_jacobian = evaluator.jacobian(initial, free_slots, (lower, upper))
    # Sensitivity counts in the coordinates the optimizer varies: with a transform, a slot whose own
    # value does not move the charges can still move them through the slots that depend on it.
    start_derivative = coordinates.derivative(internal_initial)[np.ix_(free_slots, free_slots)]
    internal_start_jacobian = start_jacobian @ start_derivative
    insensitive = [
        slot
        for column, slot in enumerate(free_slots)
        if not internal_start_jacobian[:, column].any()
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
    failed_messages: list[str] = []

    def internal_of(free_values: NDArray[np.float64]) -> NDArray[np.float64]:
        internal = internal_initial.copy()
        internal[free_slots] = free_values
        return internal

    def values_of(free_values: NDArray[np.float64]) -> NDArray[np.float64]:
        return coordinates.to_values(internal_of(free_values))

    def residuals(free_values: NDArray[np.float64]) -> NDArray[np.float64]:
        nonlocal failed
        try:
            charges = evaluator.charges(values_of(free_values))
        except EvaluationError as error:
            # Non-finite residuals make the TRF method reject the trial and shorten its step.
            failed += 1
            if len(failed_messages) < 20 and str(error) not in failed_messages:
                failed_messages.append(str(error))
            return np.full(reference.size, np.nan)
        return weight * (charges - reference)

    def jacobian(free_values: NDArray[np.float64]) -> NDArray[np.float64]:
        internal = internal_of(free_values)
        values = coordinates.to_values(internal)
        if np.array_equal(internal, internal_initial):
            charges_jacobian = start_jacobian
        else:
            charges_jacobian = evaluator.jacobian(values, free_slots, (lower, upper))
        # The chain rule from values to the internal coordinates the optimizer varies.
        derivative = coordinates.derivative(internal)[np.ix_(free_slots, free_slots)]
        return weight[:, np.newaxis] * (charges_jacobian @ derivative)

    solution = least_squares(
        residuals,
        internal_initial[free_slots],
        jac=jacobian,
        bounds=(internal_lower[free_slots], internal_upper[free_slots]),
        method="trf",
        x_scale="jac",
        max_nfev=max_evaluations,
    )
    fitted = values_of(solution.x)
    fitted.setflags(write=False)
    initial.setflags(write=False)
    final_charges = evaluator.charges(fitted)
    final_jacobian = weight[:, np.newaxis] * evaluator.jacobian(fitted, free_slots, (lower, upper))
    diagnostics = _diagnostics(
        final_charges,
        reference,
        final_jacobian,
        internal_of(solution.x),
        free_slots,
        internal_lower,
        internal_upper,
        failed_messages,
        evaluator._charge_elements(),
        tuple(evaluator.matches[slot] for slot in free_slots),
    )
    return FitResult(
        values=fitted,
        start=initial,
        free=tuple(free_slots),
        bounds=(_read_only(lower), _read_only(upper)),
        success=bool(solution.success),
        message=str(solution.message),
        evaluations=int(solution.nfev),
        failed_evaluations=failed,
        rmsd=float(np.sqrt(np.mean((final_charges - reference) ** 2))),
        diagnostics=diagnostics,
        _parameter_set=evaluator._parameter_set,
    )


def _diagnostics(
    charges: NDArray[np.float64],
    reference: NDArray[np.float64],
    jacobian: NDArray[np.float64],
    values: NDArray[np.float64],
    free: Sequence[int],
    lower: NDArray[np.float64],
    upper: NDArray[np.float64],
    failed_messages: Sequence[str],
    elements: NDArray[np.int64],
    matches: tuple[int | None, ...],
) -> Diagnostics:
    by_element = {
        periodic_table[int(atomic_number)].symbol: ElementDiagnostics(
            count=int(np.count_nonzero(elements == atomic_number)),
            rmsd=_rmsd(charges[elements == atomic_number], reference[elements == atomic_number]),
            r2=_r2(charges[elements == atomic_number], reference[elements == atomic_number]),
        )
        for atomic_number in np.unique(elements)
    }
    r2 = _r2(charges, reference)

    def at_bound(value: float, bound: float) -> bool:
        return math.isfinite(bound) and abs(value - bound) <= 1e-8 * max(1.0, abs(bound))

    at_bounds = tuple(
        slot
        for slot in free
        if at_bound(values[slot], lower[slot]) or at_bound(values[slot], upper[slot])
    )

    norms = np.linalg.norm(jacobian, axis=0)
    scaled = jacobian / np.where(norms > 0.0, norms, 1.0)
    # The right singular vectors give the weakest direction. A reduced SVD has all of them when there
    # are at least as many charges as free slots; only an underdetermined Jacobian needs the full one
    # for its null directions, and then the left vectors are small. A full SVD of a tall Jacobian
    # would build a square matrix in the number of charges.
    rows, columns = scaled.shape
    _, singular_values, right = np.linalg.svd(scaled, full_matrices=rows < columns)
    singular_values = np.concatenate((singular_values, np.zeros(len(free) - singular_values.size)))
    rank = int(np.sum(singular_values > 1e-6 * singular_values[0])) if singular_values[0] else 0
    weakest_direction = np.abs(right[-1])
    weakest = tuple(
        free[index]
        for index in np.argsort(-weakest_direction, kind="stable")
        if weakest_direction[index] >= 0.5 * weakest_direction.max()
    )
    singular_values.setflags(write=False)
    return Diagnostics(
        r2=r2,
        by_element=MappingProxyType(by_element),
        matches=matches,
        at_bounds=at_bounds,
        failed_messages=tuple(failed_messages),
        singular_values=singular_values,
        rank=rank,
        weakest=weakest,
    )


def _rmsd(charges: NDArray[np.float64], reference: NDArray[np.float64]) -> float:
    return float(np.sqrt(np.mean((charges - reference) ** 2)))


def _r2(charges: NDArray[np.float64], reference: NDArray[np.float64]) -> float:
    total = float(np.sum((reference - reference.mean()) ** 2))
    residual = float(np.sum((charges - reference) ** 2))
    return 1.0 - residual / total if total > 0.0 else math.nan


def _fixed_slots(layout: Sequence[Slot], profile: profiles.Profile | None) -> set[int]:
    if profile is None:
        return set()
    fixed = set()
    for name in profile.fixed:
        first = next((index for index, slot in enumerate(layout) if slot.name == name), None)
        if first is not None:
            fixed.add(first)
    return fixed


def _profile_bounds(
    layout: Sequence[Slot], profile: profiles.Profile | None
) -> tuple[NDArray[np.float64], NDArray[np.float64]] | None:
    if profile is None:
        return None
    lower = np.full(len(layout), -math.inf)
    upper = np.full(len(layout), math.inf)
    for index, slot in enumerate(layout):
        if slot.name in profile.bounds:
            lower[index], upper[index] = profile.bounds[slot.name]
    return lower, upper


def _start_values(
    evaluator: Evaluator,
    start: ArrayLike
    | Mapping[tuple[AtomParameterKey | BondParameterKey | None, str], float]
    | None,
    profile: profiles.Profile | None,
) -> tuple[NDArray[np.float64], set[int]]:
    """The start values, and the slots that have none."""

    if start is not None and not isinstance(start, Mapping):
        return evaluator._vector(start).copy(), set()

    values = np.array(evaluator.values)
    missing: set[int] = set()
    if profile is not None:
        groups: dict[tuple[ParameterScope, str], list[int]] = {}
        for index, slot in enumerate(evaluator.layout):
            groups.setdefault((slot.scope, slot.name), []).append(index)
        for (_, name), slots in groups.items():
            if np.all(values[slots] == 0.0):
                if name in profile.start:
                    values[slots] = profile.start[name]
                else:
                    missing.update(slots)
    if start is not None:
        positions = {(slot.key, slot.name): index for index, slot in enumerate(evaluator.layout)}
        for (key, name), value in start.items():
            if (key, name) not in positions:
                raise KeyError(f"no slot for parameter {name!r} of {key!r} in the layout")
            values[positions[(key, name)]] = value
            missing.discard(positions[(key, name)])
    return evaluator._vector(values).copy(), missing


def _read_only(values: NDArray[np.float64]) -> NDArray[np.float64]:
    copy = np.array(values)
    copy.setflags(write=False)
    return copy


def _check_start_in_bounds(
    describe: Callable[[int], str],
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
                f"slot {slot} ({describe(slot)}) = {start[slot]!r} in "
                f"[{lower[slot]!r}, {upper[slot]!r}]"
                for slot in outside
            )
        )


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
    ) -> None:
        """A transformed slot is unbounded, and is free whenever a slot it depends on is free."""

        bounded = [s for s in self._owner if math.isfinite(lower[s]) or math.isfinite(upper[s])]
        if bounded:
            raise ValueError(
                "bounds must leave these slots unbounded, because the profile transforms them and "
                "the transform sets their bounds: "
                + "; ".join(f"slot {s} ({self.describe(evaluator)(s)})" for s in bounded)
            )
        free_set = set(free)
        held = [
            s
            for s, depends_on in self._dependencies.items()
            if s not in free_set and any(d in free_set for d in depends_on)
        ]
        if held:
            raise ValueError(
                "these slots are not free, but slots their transform depends on are, so their "
                "values would change anyway; free them too, or fix the slots they depend on: "
                + "; ".join(f"slot {s} ({self.describe(evaluator)(s)})" for s in held)
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

    def describe(self, evaluator: Evaluator) -> Callable[[int], str]:
        def described(slot: int) -> str:
            text = evaluator._describe(slot)
            if slot in self._owner:
                text += f", fitted as {self._owner[slot].name}"
            return text

        return described


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
