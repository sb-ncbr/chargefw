"""Synchronous calculation facade."""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, replace
from threading import Lock
from typing import cast

from ._calculation_options import FixedAtomCharge, RequestedCalculation
from ._calculation_values import (
    CalculationCancelledError,
    CalculationResult,
    CalculationTimings,
    ChargeFWError,
    ExecutedPlan,
    Execution,
    FixedIons,
    NoExecutablePlanError,
    NumericalFailureError,
    Plan,
    Rejection,
    ResourceEstimate,
    _rejections,
)
from ._chargefw import calculation as _native_calculation
from ._chargefw import parameters as _native_parameters
from ._methods import Method, _method_catalog
from ._parameters import ParameterSet, ParameterSetCatalog, _parameter_set
from ._payloads import CalculationProgressPayload
from ._types import (
    CalculationPhase,
    ExecutionMode,
    MethodOptionValue,
    MethodScope,
    ParameterMatching,
    RequestedExecution,
)
from .core import Molecule, MoleculeCollection

__all__ = [
    "Assessment",
    "CalculationObserver",
    "CalculationProgress",
    "CalculationCancelledError",
    "CalculationResult",
    "CalculationTimings",
    "ChargeFWError",
    "ExecutedPlan",
    "FixedAtomCharge",
    "FixedIons",
    "Execution",
    "NoExecutablePlanError",
    "NumericalFailureError",
    "Plan",
    "Rejection",
    "RequestedCalculation",
    "ResourceEstimate",
    "assess",
    "calculate",
    "methods",
    "parameter_sets",
]


@dataclass(frozen=True, slots=True)
class CalculationProgress:
    """Owned snapshot of one native calculation progress event."""

    phase: CalculationPhase
    mode: ExecutionMode
    method_id: str
    target_index: int
    target_count: int
    completed_fragment_count: int
    fragment_count: int
    molecule_index: int
    conformer_index: int | None
    elapsed_seconds: float


class CalculationObserver:
    """Receives progress events and optionally requests cooperative cancellation."""

    def on_progress(self, progress: CalculationProgress) -> None:
        """Handle a progress snapshot. Implementations must be thread-safe."""

    def cancelled(self) -> bool:
        """Return whether the running calculation should stop."""

        return False


class _NativeObserverAdapter:
    __slots__ = ("_observer",)

    def __init__(self, observer: CalculationObserver) -> None:
        self._observer = observer

    def on_progress(self, progress: CalculationProgressPayload) -> None:
        self._observer.on_progress(CalculationProgress(**progress))

    def cancelled(self) -> bool:
        return self._observer.cancelled()


def _observer_adapter(observer: CalculationObserver | None) -> _NativeObserverAdapter | None:
    if observer is None:
        return None
    if not isinstance(observer, CalculationObserver):
        raise TypeError("observer must be a CalculationObserver or None")
    return _NativeObserverAdapter(observer)


for _value_type in (
    RequestedCalculation,
    Execution,
    ExecutedPlan,
    CalculationTimings,
    ChargeFWError,
    NoExecutablePlanError,
    NumericalFailureError,
    CalculationCancelledError,
    CalculationResult,
    Plan,
    Rejection,
    ResourceEstimate,
    FixedIons,
    FixedAtomCharge,
):
    _value_type.__module__ = __name__

_bundled_parameter_catalog: _native_parameters._NativeParameterCatalog | None = None
_bundled_parameter_descriptors: ParameterSetCatalog | None = None
_bundled_parameter_catalog_lock = Lock()


def _default_parameter_catalog() -> _native_parameters._NativeParameterCatalog:
    global _bundled_parameter_catalog
    if _bundled_parameter_catalog is None:
        with _bundled_parameter_catalog_lock:
            if _bundled_parameter_catalog is not None:
                return _bundled_parameter_catalog
            _bundled_parameter_catalog = _native_parameters._load_bundled_parameter_catalog()
    return _bundled_parameter_catalog


def _default_parameter_descriptors() -> ParameterSetCatalog:
    global _bundled_parameter_descriptors
    catalog = _default_parameter_catalog()
    if _bundled_parameter_descriptors is None:
        with _bundled_parameter_catalog_lock:
            if _bundled_parameter_descriptors is None:
                _bundled_parameter_descriptors = ParameterSetCatalog(
                    tuple(_parameter_set(value) for value in catalog._descriptors())
                )
    return _bundled_parameter_descriptors


parameter_sets = _default_parameter_descriptors()
methods = _method_catalog(parameter_sets)


def _as_collection(value: Molecule | MoleculeCollection | Iterable[Molecule]) -> MoleculeCollection:
    if isinstance(value, MoleculeCollection):
        return value
    if isinstance(value, Molecule):
        return MoleculeCollection((value,))
    return MoleculeCollection(value)


class Assessment:
    """Priority-ordered runnable plans and rejected alternatives for immutable molecules."""

    __slots__ = ("_native", "_requested", "_plans", "_rejections", "_seconds")

    def __init__(
        self,
        native: _native_calculation._NativeAssessment,
        molecules: MoleculeCollection,
        requested: RequestedCalculation,
    ) -> None:
        payload = native.report()
        self._native = native
        self._requested = requested
        self._plans = tuple(
            Plan(plan, molecules, methods, parameter_sets, requested) for plan in native.plans()
        )
        self._rejections = _rejections(payload["rejections"], methods, parameter_sets)
        self._seconds = payload["applicability_seconds"]

    @property
    def plans(self) -> tuple[Plan, ...]:
        return self._plans

    @property
    def rejections(self) -> tuple[Rejection, ...]:
        return self._rejections

    @property
    def default_plan(self) -> Plan | None:
        return None if not self._plans else self._plans[0]

    @property
    def seconds(self) -> float:
        return self._seconds


def assess(
    molecules: Molecule | MoleculeCollection | Iterable[Molecule],
    *,
    method: str | Method | None = None,
    parameter_set: str | ParameterSet | None = None,
    method_scope: MethodScope = "recommended",
    options: Mapping[str, MethodOptionValue] | None = None,
    options_by_method: Mapping[str, Mapping[str, MethodOptionValue]] | None = None,
    parameter_matching: ParameterMatching = "strict",
    execution: RequestedExecution = "auto",
    radius: float | None = None,
    cutoff_threshold: int | None = 20_000,
    cover_threshold: int | None = 80_000,
    threads: int = 0,
    fixed_ions: Sequence[str] | None = None,
) -> Assessment:
    """Inspect applicability and return reusable executable calculation plans."""

    requested = RequestedCalculation(
        method=method,
        parameter_set=parameter_set,
        method_scope=method_scope,
        options=options,
        options_by_method=options_by_method,
        parameter_matching=parameter_matching,
        execution=execution,
        radius=radius,
        cutoff_threshold=cutoff_threshold,
        cover_threshold=cover_threshold,
        threads=threads,
        fixed_ions=fixed_ions,
    )
    collection = _as_collection(molecules)
    native = _native_calculation._make_assessment(
        collection._native_molecules,
        tuple(getattr(molecule, "_native_input_metadata", None) for molecule in collection),
        tuple(
            (
                molecule.source_name,
                molecule.record_index,
                molecule.record_id,
            )
            for molecule in collection
        ),
        tuple(molecule._caller_atom_ids for molecule in collection),
        collection.name,
        _default_parameter_catalog(),
        requested.method,
        requested.parameter_set,
        requested.method_scope,
        {
            method_id: dict(overrides)
            for method_id, overrides in requested.options_by_method.items()
        },
        requested._permissive_types,
        requested.execution,
        requested.radius,
        requested.cutoff_threshold,
        requested.cover_threshold,
        requested.threads,
        requested.fixed_ions,
    )
    return Assessment(native, collection, requested)


class _OmittedPlan:
    def __repr__(self) -> str:
        return "<automatic>"


_OMITTED_PLAN = _OmittedPlan()


class _Default:
    """A default argument value that records that the caller did not pass the argument."""

    __slots__ = ("value",)

    def __init__(self, value: object) -> None:
        self.value = value

    def __repr__(self) -> str:
        return repr(self.value)


def _passed[T](value: T) -> bool:
    return not isinstance(value, _Default)


def _value[T](value: T) -> T:
    return cast(T, value.value) if isinstance(value, _Default) else value


def calculate(
    molecules: Molecule | MoleculeCollection | Iterable[Molecule],
    plan: Plan = cast(Plan, _OMITTED_PLAN),
    *,
    method: str | Method | None = None,
    parameter_set: str | ParameterSet | None = None,
    method_scope: MethodScope = cast(MethodScope, _Default("recommended")),
    options: Mapping[str, MethodOptionValue] | None = None,
    options_by_method: Mapping[str, Mapping[str, MethodOptionValue]] | None = None,
    parameter_matching: ParameterMatching = cast(ParameterMatching, _Default("strict")),
    execution: RequestedExecution = cast(RequestedExecution, _Default("auto")),
    radius: float | None = None,
    cutoff_threshold: int | None = cast(int, _Default(20_000)),
    cover_threshold: int | None = cast(int, _Default(80_000)),
    threads: int | None = None,
    observer: CalculationObserver | None = None,
    fixed_ions: Sequence[str] | None = None,
) -> CalculationResult:
    """Calculate molecules directly, or execute an explicitly assessed plan.

    Direct calculation assesses the input and uses the default plan. Without an explicit method,
    selection uses the deterministic catalog order among applicable methods.
    """

    collection = _as_collection(molecules)
    native_observer = _observer_adapter(observer)
    if isinstance(plan, _OmittedPlan):
        assessment = assess(
            collection,
            method=method,
            parameter_set=parameter_set,
            method_scope=_value(method_scope),
            options=options,
            options_by_method=options_by_method,
            parameter_matching=_value(parameter_matching),
            execution=_value(execution),
            radius=radius,
            cutoff_threshold=_value(cutoff_threshold),
            cover_threshold=_value(cover_threshold),
            threads=0 if threads is None else threads,
            fixed_ions=fixed_ions,
        )
        if assessment.default_plan is None:
            result = CalculationResult(
                assessment._native.calculate_default(native_observer),
                collection,
                assessment._requested,
                methods,
                parameter_sets,
                assessment.rejections,
            )
            result._raise_for_status()
            return result
        return calculate(collection, assessment.default_plan, threads=threads, observer=observer)

    if plan is None or not isinstance(plan, Plan):
        raise TypeError("plan must be a Plan; omit it to use automatic assessment")
    if any(
        (
            method is not None,
            parameter_set is not None,
            _passed(method_scope),
            options is not None,
            options_by_method is not None,
            _passed(parameter_matching),
            _passed(execution),
            radius is not None,
            _passed(cutoff_threshold),
            _passed(cover_threshold),
            fixed_ions is not None,
        )
    ):
        raise TypeError("selection arguments cannot be combined with an assessed plan")
    if not plan._matches(collection):
        raise ValueError("plan was assessed for a different molecule collection")
    normalized_threads = None if threads is None else RequestedCalculation(threads=threads).threads
    requested = (
        plan._requested
        if normalized_threads is None
        else replace(plan._requested, threads=normalized_threads)
    )
    result = CalculationResult(
        plan._calculate(normalized_threads, native_observer),
        collection,
        requested,
        methods,
        parameter_sets,
    )
    result._raise_for_status()
    return result


Assessment.__module__ = __name__
