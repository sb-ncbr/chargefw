"""Fitting parameter values to reference charges in the least-squares sense."""

from __future__ import annotations

import math
import time
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import Literal, SupportsIndex

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .._arrays import immutable_array
from .._parameters import AtomParameterKey, BondParameterKey, ParameterSet
from .._types import ParameterScope
from . import profiles
from ._compression import CompressionError, compress, determined_step
from ._coordinates import _Coordinates
from ._diagnostics import Diagnostics, _diagnostics, _spectrum_of_gram, _spectrum_of_jacobian
from ._evaluator import EvaluationError, Evaluator, Slot, _slot_values
from ._observer import FitObserver, FitProgress
from ._template import _with_values


@dataclass(frozen=True, slots=True, eq=False)
class FitResult:
    """The outcome of ``fit()``.

    ``values`` and ``start`` hold every slot in the order of the evaluator's layout; ``free`` lists
    the slots the fit changed, and ``bounds`` the lower and upper values it used per slot; a slot the
    profile transforms has no bounds of its own there, as its bound applies to the transformed
    quantity. ``rmsd`` compares the charges at ``values`` with the references, without weights,
    leaving out charges without a reference. ``failed_evaluations`` counts trial values whose
    calculation failed and were therefore rejected.
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

        return _with_values(self._parameter_set, self.values, id, name, publication, notes)


class _Cancelled(Exception):
    """Raised inside an evaluation when the observer asks the fit to stop."""


def fit(
    evaluator: Evaluator,
    references: Sequence[ArrayLike] | ArrayLike,
    *,
    weights: Sequence[ArrayLike] | ArrayLike | None = None,
    start: ArrayLike
    | Mapping[tuple[AtomParameterKey | BondParameterKey | None, str], float]
    | ParameterSet
    | None = None,
    free: Sequence[SupportsIndex] | None = None,
    bounds: tuple[ArrayLike, ArrayLike] | None = None,
    profile: profiles.Profile | Literal["default"] | None = "default",
    max_evaluations: int | None = None,
    observer: FitObserver | None = None,
) -> FitResult:
    """Fit parameter values so that the charges match reference charges in the least-squares sense.

    ``references`` holds the reference charges of every assignment, in assignment order: one array
    per assignment, or one array for all. A reference of ``nan`` means the charge has none: it does
    not enter the fit, whatever its weight, and the RMSD and the reference-based error statistics of
    the diagnostics leave it out.
    ``weights`` has the same shape, or is omitted for equal weights. ``start`` is a vector in layout
    order; a mapping from (key, parameter name) to the values that differ from the evaluator's
    values, the key being ``None`` for common parameters; or a parameter set, such as an earlier
    fit or a checkpoint, whose values are used for every slot with the same key and parameter name,
    for any method and layout, a bond in either endpoint order; slots it lacks start as without
    ``start``. This copies values; it does not resume an earlier fit, and the profile's start of the
    linear parameters applies as with any other start. ``free`` lists the slots to fit; the others
    keep their start values. ``bounds`` gives lower and upper values per slot.

    ``observer`` receives a ``FitProgress`` after every trial point the optimizer evaluates, failed
    ones included; not during the checks and the start before the first trial point, the derivatives
    of a method without exact ones, or the final diagnostics. When it asks the fit to stop, the fit
    calculates the diagnostics at the best values so far, which costs one more evaluation (and, for a
    method without exact derivatives, its forward differences), and returns them with ``success``
    false. If that calculation fails, the fit raises ``EvaluationError``; the observer has received
    the best values before.

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

    The fit takes SciPy's trust-region reflective steps. For a method with exact derivatives (EEM,
    SQE, SQE+q0, SQE+qp), each step uses them in the form of the normal equations, whose size
    depends on the free slots and not on the number of charges. The fit then does not try to improve
    the charges through combinations of free slots that the data barely determines, which the
    diagnostics usually report through ``rank``; plain least squares would follow them, usually to
    extreme values. Other methods use the forward differences of ``Evaluator.jacobian()``.

    Before fitting, the start is checked: it must lie within the bounds of every free slot; the
    charges must be calculable; and every free slot must change a charge whose weight is not zero. A
    start on a bound, as from an earlier fit that ended there, or beyond it by no more than
    ``1e-10 * max(1, |bound|)``, is accepted; SciPy moves it just inside the bound.
    Calculations that fail at trial values are rejected, and the optimizer tries a shorter step. A
    finite-difference step that fails is tried the other way; if that fails too, the fit stops with
    ``EvaluationError``. So does a fit whose normal equations cannot represent its least-squares
    model in floating point. Needs SciPy.
    """

    started = time.perf_counter()
    try:
        from scipy.optimize import least_squares  # type: ignore[import-untyped, unused-ignore]
    except ImportError as error:
        raise ImportError(
            "fit() needs SciPy; install it, for example with pip install chargefw[parameterization]"
        ) from error

    chosen = (
        profiles.PROFILES.get(evaluator._parameter_set.method) if profile == "default" else profile
    )
    layout = evaluator.layout
    lower, upper = evaluator._bounds(
        bounds if bounds is not None else _profile_bounds(layout, chosen)
    )
    coordinates = _Coordinates(chosen.transforms if chosen is not None else (), layout)
    # The fit works in internal coordinates: the values, except that transformed slots hold the
    # transformed quantity, which carries the transform's bounds.
    internal_lower, internal_upper = coordinates.bounds(lower, upper)
    # Slots the profile holds at the lower bound of what the fit varies, each with the slots that
    # take up the difference; only by default.
    compensation = _terminal_slots(evaluator, chosen, internal_lower, upper) if free is None else {}
    held = set(compensation)
    if free is not None:
        free_slots = evaluator._slot_indices(free)
    else:
        fixed = _fixed_slots(layout, chosen) | held
        free_slots = [slot for slot in range(len(layout)) if slot not in fixed]
    if not free_slots:
        raise ValueError("fit needs at least one free slot")
    if len(set(free_slots)) != len(free_slots):
        raise ValueError("free slots must not repeat")
    coordinates.check(evaluator, free_slots, lower, upper, held)
    initial, missing = _start_values(evaluator, start, chosen)
    missing_free = [slot for slot in free_slots if slot in missing]
    if missing_free:
        raise ValueError(
            "these free slots have no start value, because their parameter is not given in the "
            "parameter set and has no default; give them in start: "
            + "; ".join(f"slot {slot} ({evaluator._describe(slot)})" for slot in missing_free)
        )
    internal_initial = coordinates.to_internal(initial)
    if held:
        # Each held value goes onto its limit, and its compensating slots change by the difference,
        # so the start keeps its charges. The result reports this start.
        slots = sorted(held)
        before = initial[slots].copy()
        internal_initial[slots] = internal_lower[slots]
        initial[slots] = coordinates.to_values(internal_initial)[slots]
        for slot, old in zip(slots, before, strict=True):
            for other in compensation[slot]:
                initial[other] += old - initial[slot]
        internal_initial = coordinates.to_internal(initial)
        internal_initial[slots] = internal_lower[slots]
    _check_start_in_bounds(
        coordinates.describe(evaluator),
        internal_initial,
        free_slots,
        internal_lower,
        internal_upper,
    )
    # A start on a bound, or beyond it by rounding, goes onto the bound, which every calculation from
    # here on uses; SciPy then moves it just inside. The result reports the start as given.
    internal_initial[free_slots] = np.clip(
        internal_initial[free_slots], internal_lower[free_slots], internal_upper[free_slots]
    )
    feasible = coordinates.to_values(internal_initial)

    try:
        start_charges = evaluator.charges(feasible)
    except EvaluationError as error:
        raise EvaluationError(f"the charges cannot be calculated at the start: {error}") from error
    assert evaluator._assignment_sizes is not None
    sizes = evaluator._assignment_sizes
    reference = _flat(references, "references", sizes, allow_missing=True)
    weight = np.ones(start_charges.size) if weights is None else _flat(weights, "weights", sizes)
    # A charge without a reference gets weight zero; `target` holds a placeholder there, which the
    # weight removes from every sum the fit forms. `reference` keeps the nan for the diagnostics.
    present = ~np.isnan(reference)
    weight = np.where(present, weight, 0.0)
    target = np.where(present, reference, 0.0)

    # By default, the free slots of the profile's linear parameters start from their best values for
    # the start of the others; the result reports this start.
    linear = (
        [
            slot
            for slot in free_slots
            if layout[slot].name in chosen.linear
            and coordinates.plain(slot)
            and np.isneginf(lower[slot])
            and np.isposinf(upper[slot])
        ]
        if free is None and chosen is not None
        else []
    )
    if linear:
        try:
            feasible[linear] += _linear_step(
                evaluator, feasible, linear, target, weight, (lower, upper)
            )
        except CompressionError as error:
            raise EvaluationError(
                f"the best start of the linear parameters cannot be calculated: {error}"
            ) from error
        initial[linear] = feasible[linear]
        internal_initial[linear] = feasible[linear]

    exact = evaluator._exact
    free_index = np.ix_(free_slots, free_slots)
    # The values that change with the free internal coordinates: the free slots, and the slots the
    # profile holds in internal coordinates, whose values follow the free slots they depend on. The
    # chain rule maps derivatives in these values onto the free internal coordinates.
    moving = free_slots + sorted(held)
    moving_index = np.ix_(moving, moving)
    chain_index = np.ix_(moving, free_slots)
    # Sensitivity counts in the coordinates the optimizer varies: with a transform, a slot whose own
    # value does not move the charges can still move them through the slots that depend on it.
    # Charges with weight zero do not count.
    start_derivative = coordinates.derivative(internal_initial)[chain_index]
    if exact:
        _, _, start_gram, _ = evaluator._normal_equations(feasible, target, weight)
        sensitivity = np.diag(start_derivative.T @ start_gram[moving_index] @ start_derivative)
        sensitive = sensitivity != 0.0
    else:
        start_jacobian = evaluator.jacobian(feasible, moving, (lower, upper))
        internal_start_jacobian = (weight[:, np.newaxis] * start_jacobian) @ start_derivative
        sensitive = internal_start_jacobian.any(axis=0)
    insensitive = [slot for column, slot in enumerate(free_slots) if not sensitive[column]]
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
    count = 0
    best_objective = math.inf
    best_rmsd = math.nan
    best_free: NDArray[np.float64] | None = None

    def internal_of(free_values: NDArray[np.float64]) -> NDArray[np.float64]:
        internal = internal_initial.copy()
        internal[free_slots] = free_values
        return internal

    def values_of(free_values: NDArray[np.float64]) -> NDArray[np.float64]:
        return coordinates.to_values(internal_of(free_values))

    def record(error: EvaluationError) -> None:
        nonlocal failed
        failed += 1
        if len(failed_messages) < 20 and str(error) not in failed_messages:
            failed_messages.append(str(error))

    def observe(
        free_values: NDArray[np.float64], objective: float, charges: NDArray[np.float64] | None
    ) -> None:
        """Keeps the best point, tells the observer, and stops the fit when it asks to."""

        nonlocal count, best_objective, best_rmsd, best_free
        count += 1
        rmsd = math.nan if charges is None else _rmsd(charges, reference)
        if math.isfinite(objective) and objective < best_objective:
            best_objective, best_rmsd, best_free = objective, rmsd, free_values.copy()
        if observer is None:
            return
        best_values = immutable_array(values_of(best_free) if best_free is not None else initial)
        observer.on_evaluation(
            FitProgress(
                evaluation=count,
                failed_evaluations=failed,
                objective=objective,
                rmsd=rmsd,
                best_objective=best_objective,
                best_rmsd=best_rmsd,
                best_values=best_values,
                elapsed_seconds=time.perf_counter() - started,
                _parameter_set=evaluator._parameter_set,
            )
        )
        if observer.cancelled():
            raise _Cancelled

    if exact:
        # The least-squares model in the size of the free slots, from the normal equations; see
        # _compression. Cached, so that the residuals and the matrix of one point share one call.
        size = len(free_slots)
        cache: dict[bytes, tuple[NDArray[np.float64], NDArray[np.float64]]] = {}

        def compressed(
            free_values: NDArray[np.float64],
        ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
            key = free_values.tobytes()
            if key not in cache:
                internal = internal_of(free_values)
                try:
                    charges, total, gram, gradient = evaluator._normal_equations(
                        coordinates.to_values(internal), target, weight
                    )
                except EvaluationError as error:
                    # Non-finite residuals make the TRF method reject the trial and shorten its
                    # step.
                    record(error)
                    observe(free_values, math.nan, None)
                    represented = (np.full(size + 1, np.nan), np.zeros((size + 1, size)))
                else:
                    observe(free_values, total, charges)
                    derivative = coordinates.derivative(internal)[chain_index]
                    try:
                        represented = compress(
                            total,
                            derivative.T @ gram[moving_index] @ derivative,
                            derivative.T @ gradient[moving],
                        )
                    except CompressionError as error:
                        raise EvaluationError(
                            f"the fit cannot continue with exact derivatives: {error}"
                        ) from error
                cache.clear()
                cache[key] = represented
            return cache[key]

        def residuals(free_values: NDArray[np.float64]) -> NDArray[np.float64]:
            return compressed(free_values)[0]

        def jacobian(free_values: NDArray[np.float64]) -> NDArray[np.float64]:
            return compressed(free_values)[1]

    else:

        def residuals(free_values: NDArray[np.float64]) -> NDArray[np.float64]:
            try:
                charges = evaluator.charges(values_of(free_values))
            except EvaluationError as error:
                # Non-finite residuals make the TRF method reject the trial and shorten its step.
                record(error)
                observe(free_values, math.nan, None)
                return np.full(reference.size, np.nan)
            residual = weight * (charges - target)
            observe(free_values, float(residual @ residual), charges)
            return residual

        def jacobian(free_values: NDArray[np.float64]) -> NDArray[np.float64]:
            internal = internal_of(free_values)
            values = coordinates.to_values(internal)
            if np.array_equal(internal, internal_initial):
                charges_jacobian = start_jacobian
            else:
                charges_jacobian = evaluator.jacobian(values, moving, (lower, upper))
            # The chain rule from values to the internal coordinates the optimizer varies.
            derivative = coordinates.derivative(internal)[chain_index]
            return weight[:, np.newaxis] * (charges_jacobian @ derivative)

    try:
        solution = least_squares(
            residuals,
            internal_initial[free_slots],
            jac=jacobian,
            bounds=(internal_lower[free_slots], internal_upper[free_slots]),
            method="trf",
            x_scale="jac",
            loss="linear",
            tr_solver="exact",
            max_nfev=max_evaluations,
        )
        solved, success, message, evaluations = (
            solution.x,
            bool(solution.success),
            str(solution.message),
            int(solution.nfev),
        )
    except _Cancelled:
        if best_free is None:
            raise EvaluationError("the fit was cancelled before any evaluation succeeded") from None
        solved, success, message, evaluations = (
            best_free,
            False,
            "cancelled by the observer; the values are the best so far",
            count,
        )
    fitted = immutable_array(values_of(solved))
    initial = immutable_array(initial)
    # The diagnostics use the weighted derivatives in the parameter values, not the internal
    # coordinates.
    if exact:
        final_charges, _, final_gram, _ = evaluator._normal_equations(fitted, target, weight)
        spectrum = _spectrum_of_gram(final_gram[free_index])
    else:
        final_charges = evaluator.charges(fitted)
        spectrum = _spectrum_of_jacobian(
            weight[:, np.newaxis] * evaluator.jacobian(fitted, free_slots, (lower, upper))
        )
    diagnostics = _diagnostics(
        final_charges,
        reference,
        spectrum,
        internal_of(solved),
        free_slots,
        internal_lower,
        internal_upper,
        failed_messages,
        evaluator._charge_elements(),
        evaluator._charge_atom_entries(),
        evaluator._parameter_set.atom_keys,
        tuple(evaluator.matches[slot] for slot in free_slots),
    )
    return FitResult(
        values=fitted,
        start=initial,
        free=tuple(free_slots),
        bounds=(immutable_array(lower), immutable_array(upper)),
        success=success,
        message=message,
        evaluations=evaluations,
        failed_evaluations=failed,
        rmsd=_rmsd(final_charges, reference),
        diagnostics=diagnostics,
        _parameter_set=evaluator._parameter_set,
    )


def _terminal_slots(
    evaluator: Evaluator,
    profile: profiles.Profile | None,
    internal_lower: NDArray[np.float64],
    upper: NDArray[np.float64],
) -> dict[int, tuple[int, ...]]:
    """The slots the profile's ``terminal`` holds, each with the bond slots that take up the
    difference: those of its type's bonds with the named parameter. A slot is held only when the
    lower bound of what the fit varies for it is finite and none of those bond slots has an upper
    bound, so that the difference can always be taken up."""

    if profile is None:
        return {}
    bond_slots: dict[tuple[int, str], int] = {}
    for index, slot in enumerate(evaluator.layout):
        if slot.scope == "bond":
            bond_slots[(evaluator._entry_of[index], slot.name)] = index
    result = {}
    for index, slot in enumerate(evaluator.layout):
        bonds = evaluator._terminal[index]
        if slot.name not in profile.terminal or bonds is None:
            continue
        name = profile.terminal[slot.name]
        others = tuple(bond_slots[(bond, name)] for bond in bonds if (bond, name) in bond_slots)
        if (
            len(others) == len(bonds)
            and math.isfinite(internal_lower[index])
            and not any(math.isfinite(upper[other]) for other in others)
        ):
            result[index] = others
    return result


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
    | ParameterSet
    | None,
    profile: profiles.Profile | None,
) -> tuple[NDArray[np.float64], set[int]]:
    """The start values, and the slots that have none."""

    if isinstance(start, ParameterSet):
        # The values of the slots the set shares with the layout; the set's other entries do not
        # matter, and the layout's other slots start as without a start.
        layout_keys = {_unordered(slot.key, slot.name) for slot in evaluator.layout}
        start = {
            key: value
            for key, value in _slot_values(start).items()
            if _unordered(*key) in layout_keys
        }
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
        positions = {
            _unordered(slot.key, slot.name): index for index, slot in enumerate(evaluator.layout)
        }
        for (key, name), value in start.items():
            position = positions.get(_unordered(key, name))
            if position is None:
                raise KeyError(f"no slot for parameter {name!r} of {key!r} in the layout")
            values[position] = value
            missing.discard(position)
    return evaluator._vector(values).copy(), missing


def _unordered(
    key: AtomParameterKey | BondParameterKey | None, name: str
) -> tuple[AtomParameterKey | BondParameterKey | None, str]:
    """A slot's key and parameter name, with a bond's endpoints in a fixed order: a bond key means
    the same in either order, and parameter files keep the order they were written in."""

    if isinstance(key, BondParameterKey):
        first, second = sorted(
            key.atoms, key=lambda atom: (atom.element, atom.classifier, atom.type)
        )
        return (
            BondParameterKey(atoms=(first, second), classifier=key.classifier, type=key.type),
            name,
        )
    return (key, name)


def _linear_step(
    evaluator: Evaluator,
    values: NDArray[np.float64],
    slots: list[int],
    reference: NDArray[np.float64],
    weight: NDArray[np.float64],
    bounds: tuple[NDArray[np.float64], NDArray[np.float64]],
) -> NDArray[np.float64]:
    """The change of `slots` that minimizes the weighted sum of squares when the charges are affine
    in them, in the directions the data determines."""

    if evaluator._exact:
        _, _, gram, gradient = evaluator._normal_equations(values, reference, weight)
        return determined_step(gram[np.ix_(slots, slots)], gradient[slots])
    jacobian = weight[:, np.newaxis] * evaluator.jacobian(values, slots, bounds)
    residuals = weight * (evaluator.charges(values) - reference)
    return determined_step(jacobian.T @ jacobian, jacobian.T @ residuals)


def _check_start_in_bounds(
    describe: Callable[[int], str],
    start: NDArray[np.float64],
    free: Sequence[int],
    lower: NDArray[np.float64],
    upper: NDArray[np.float64],
) -> None:
    # A value an earlier fit left on a bound can come back across it by rounding, for example after
    # a transform; beyond that, the start is outside.
    def tolerance(bound: float) -> float:
        return 1e-10 * max(1.0, abs(bound)) if math.isfinite(bound) else 0.0

    outside = [
        slot
        for slot in free
        if not lower[slot] - tolerance(lower[slot])
        <= start[slot]
        <= upper[slot] + tolerance(upper[slot])
    ]
    if outside:
        raise ValueError(
            "the start must lie within the bounds of every free slot; change the start or the "
            "bounds of: "
            + "; ".join(
                f"slot {slot} ({describe(slot)}) = {start[slot]!r} in "
                f"[{lower[slot]!r}, {upper[slot]!r}]"
                for slot in outside
            )
        )


def _rmsd(charges: NDArray[np.float64], reference: NDArray[np.float64]) -> float:
    """The RMSD over the charges that have a reference; nan when none has."""

    present = ~np.isnan(reference)
    return (
        float(np.sqrt(np.mean((charges[present] - reference[present]) ** 2)))
        if present.any()
        else math.nan
    )


def _flat(
    value: Sequence[ArrayLike] | ArrayLike,
    name: str,
    sizes: Sequence[int],
    *,
    allow_missing: bool = False,
) -> NDArray[np.float64]:
    """One value per calculated charge, from one array per assignment or one array for all; with
    ``allow_missing``, nan marks a missing value."""

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
    if allow_missing:
        if np.any(np.isinf(flat)):
            raise ValueError(f"{name} must be finite, or nan where there is none")
        if np.all(np.isnan(flat)):
            raise ValueError(f"{name} has no value at all")
    elif not np.all(np.isfinite(flat)):
        raise ValueError(f"{name} must be finite")
    return flat
