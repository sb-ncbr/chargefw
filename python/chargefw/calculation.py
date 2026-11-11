"""Synchronous calculation facade."""

from __future__ import annotations

import os
import warnings
from collections.abc import Iterable, Mapping, Sequence
from threading import Lock

from ._calculation_options import FixedAtomCharge, RequestedCalculation
from ._calculation_values import (
    CalculationCancelledError,
    CalculationObserver,
    CalculationProgress,
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
    _observer_adapter,
    _outcome,
    _rejections,
)
from ._chargefw import calculation as _native_calculation
from ._chargefw import parameters as _native_parameters
from ._methods import Method, _method_catalog
from ._parameters import ParameterSet, ParameterSetCatalog, _parameter_set
from ._public import set_public_module
from ._types import (
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


set_public_module(
    __name__,
    CalculationProgress,
    CalculationObserver,
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


def calculate(
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
    observer: CalculationObserver | None = None,
    fixed_ions: Sequence[str] | None = None,
) -> CalculationResult:
    """Assess the molecules and calculate the default plan.

    Without an explicit method, selection uses the deterministic catalog order among applicable
    methods. To choose another plan, use ``assess()`` and ``Plan.calculate()``.
    """

    collection = _as_collection(molecules)
    native_observer = _observer_adapter(observer)
    assessment = assess(
        collection,
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
    plan = assessment.default_plan
    if plan is None:
        result = CalculationResult(
            assessment._native.calculate_default(native_observer),
            collection,
            assessment._requested,
            methods,
            assessment._parameter_sets,
            assessment.rejections,
        )
    else:
        result = plan._execute(None, native_observer)

    # Raised here rather than in a helper, so that the traceback ends at this function.
    error, cautions = _outcome(result)
    if error is not None:
        raise error
    for message, category in cautions:
        warnings.warn(message, category, stacklevel=2)
    return result


set_public_module(__name__, Assessment)
