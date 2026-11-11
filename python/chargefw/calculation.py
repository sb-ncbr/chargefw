"""Synchronous calculation facade."""

from __future__ import annotations

import os
import warnings
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
    ChargeIssue,
    ExecutedPlan,
    Execution,
    FixedIons,
    InputIssue,
    NoExecutablePlanError,
    NumericalFailureError,
    Plan,
    Rejection,
    ResourceEstimate,
    _input_issues,
    _rejections,
)
from ._chargefw import calculation as _native_calculation
from ._chargefw import parameters as _native_parameters
from ._methods import Method, _method_catalog
from ._parameters import ParameterSet, ParameterSetCatalog, _parameter_set
from ._payloads import CalculationProgressPayload
from ._public import set_public_module
from ._types import (
    CalculationPhase,
    ExecutionMode,
    MethodOptionValue,
    MethodScope,
    ParameterMatching,
    RequestedExecution,
)
from ._warnings import (
    CHARGE_WARNING_KINDS,
    INPUT_WARNING_KINDS,
    ChargeFWWarning,
    ChargeWarning,
    InputWarning,
    issue_warnings,
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
    "InputIssue",
    "ChargeIssue",
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
    "read_parameter_set",
    "write_parameter_set",
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


set_public_module(
    __name__,
    RequestedCalculation,
    Execution,
    ExecutedPlan,
    CalculationTimings,
    ChargeFWError,
    NoExecutablePlanError,
    NumericalFailureError,
    CalculationCancelledError,
    CalculationResult,
    InputIssue,
    ChargeIssue,
    Plan,
    Rejection,
    ResourceEstimate,
    FixedIons,
    FixedAtomCharge,
)

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


def _check_custom_parameter_id(id: str, source: str | None = None) -> None:
    if id in parameter_sets:
        where = f" in {source!r}" if source else ""
        raise ValueError(f"parameter set {id!r}{where} has the ID of a bundled parameter set")


def _parameter_set_from_catalog(
    native: _native_parameters._NativeParameterCatalog,
) -> ParameterSet:
    (descriptor,) = native._descriptors()
    return _parameter_set(descriptor, native)


def _native_catalog_for(parameter_set: ParameterSet) -> _native_parameters._NativeParameterCatalog:
    return parameter_set._native or _default_parameter_catalog()


def read_parameter_set(path: str | os.PathLike[str]) -> ParameterSet:
    """Read a parameter set from a JSON file in the ChargeFW parameter-set format.

    Pass the result as ``parameter_set`` to ``assess()`` or ``calculate()``, with its method. Its ID
    must differ from the IDs of the bundled parameter sets.
    """

    native = _native_parameters._load_parameter_set_catalog(os.fspath(path))
    result = _parameter_set_from_catalog(native)
    _check_custom_parameter_id(result.id, os.fspath(path))
    return result


def write_parameter_set(path: str | os.PathLike[str], parameter_set: ParameterSet) -> None:
    """Write a parameter set to a JSON file in the ChargeFW parameter-set format.

    Works for bundled sets and for sets read from files or made by fitting. The file is created or
    replaced.
    """

    if not isinstance(parameter_set, ParameterSet):
        raise TypeError("parameter_set must be a ParameterSet")

    native = _native_catalog_for(parameter_set)
    text = native._json(parameter_set.id)
    with open(path, "w", encoding="utf-8") as output:
        output.write(text)


def _parameter_catalogs(
    parameter_set: str | ParameterSet | None,
) -> tuple[_native_parameters._NativeParameterCatalog, ParameterSetCatalog]:
    """The native catalog for an assessment and the catalog that resolves its results.

    A set that is not bundled is assessed alone, and its results resolve it next to the bundled
    sets; read_parameter_set() guarantees that its ID does not collide with theirs.
    """

    if isinstance(parameter_set, ParameterSet) and parameter_set._native is not None:
        return parameter_set._native, ParameterSetCatalog((*parameter_sets.values(), parameter_set))
    return _default_parameter_catalog(), parameter_sets


def _as_collection(value: Molecule | MoleculeCollection | Iterable[Molecule]) -> MoleculeCollection:
    if isinstance(value, MoleculeCollection):
        return value
    if isinstance(value, Molecule):
        return MoleculeCollection((value,))
    return MoleculeCollection(value)


class Assessment:
    """Priority-ordered runnable plans and rejected alternatives for immutable molecules."""

    __slots__ = (
        "_native",
        "_requested",
        "_parameter_sets",
        "_plans",
        "_rejections",
        "_input_warnings",
        "_seconds",
    )

    def __init__(
        self,
        native: _native_calculation._NativeAssessment,
        molecules: MoleculeCollection,
        requested: RequestedCalculation,
        parameter_sets: ParameterSetCatalog,
    ) -> None:
        payload = native.report()
        self._native = native
        self._requested = requested
        self._parameter_sets = parameter_sets
        self._plans = tuple(
            Plan(plan, molecules, methods, parameter_sets, requested) for plan in native.plans()
        )
        self._rejections = _rejections(payload["rejections"], methods, parameter_sets)
        self._input_warnings = _input_issues(payload["input_warnings"])
        self._seconds = payload["applicability_seconds"]

    @property
    def plans(self) -> tuple[Plan, ...]:
        return self._plans

    @property
    def rejections(self) -> tuple[Rejection, ...]:
        return self._rejections

    @property
    def input_warnings(self) -> tuple[InputIssue, ...]:
        """Every input issue of the molecules, whichever plan is executed."""

        return self._input_warnings

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
    cutoff_threshold: int | None = _native_calculation._default_cutoff_atom_threshold,
    cover_threshold: int | None = _native_calculation._default_cover_atom_threshold,
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
    native_catalog, catalog = _parameter_catalogs(parameter_set)
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
        native_catalog,
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
    return Assessment(native, collection, requested, catalog)


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
    cutoff_threshold: int | None = cast(
        int, _Default(_native_calculation._default_cutoff_atom_threshold)
    ),
    cover_threshold: int | None = cast(
        int, _Default(_native_calculation._default_cover_atom_threshold)
    ),
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
                assessment._parameter_sets,
                assessment.rejections,
            )
        else:
            result = _calculate_plan(collection, assessment.default_plan, threads, native_observer)
    else:
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
        result = _calculate_plan(collection, plan, threads, native_observer)

    # Raised here rather than in a helper, so that the traceback ends at this function.
    error, cautions = _outcome(result)
    if error is not None:
        raise error
    for message, category in cautions:
        warnings.warn(message, category, stacklevel=2)
    return result


def _calculate_plan(
    collection: MoleculeCollection,
    plan: Plan,
    threads: int | None,
    native_observer: _NativeObserverAdapter | None,
) -> CalculationResult:
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
        plan._parameter_sets,
    )
    return result


def _outcome(
    result: CalculationResult,
) -> tuple[ChargeFWError | None, list[tuple[str, type[ChargeFWWarning]]]]:
    """The error of a failed calculation, or the warnings about issues that can make charges wrong.

    A failure has no warnings: its error carries the issues in its result, also when warnings are
    turned into errors.
    """

    error = result._status_error()
    if error is not None:
        return error, []
    molecules = result.molecules

    def describe(index: int) -> str:
        name = molecules[index].name
        return f"molecule {index + 1}" + (f" ({name!r})" if name else "")

    return None, [
        *issue_warnings(
            ((issue.molecule_index, issue) for issue in result.input_warnings),
            kinds=INPUT_WARNING_KINDS,
            category=InputWarning,
            molecule_count=len(molecules),
            describe=describe,
            attribute="result.input_warnings",
        ),
        *issue_warnings(
            ((issue.molecule_index, issue) for issue in result.charge_warnings),
            kinds=CHARGE_WARNING_KINDS,
            category=ChargeWarning,
            molecule_count=len(molecules),
            describe=describe,
            attribute="result.charge_warnings",
        ),
    ]


set_public_module(__name__, Assessment)
