"""Private calculation values: plans, results, issues, and observers, with payload conversion."""

from __future__ import annotations

import warnings
from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field, replace
from types import MappingProxyType

import numpy as np
from numpy.typing import NDArray

from ._arrays import as_integer, index_array
from ._calculation_options import FixedAtomCharge, RequestedCalculation
from ._chargefw import calculation as _native_calculation
from ._charges import ChargeAssignment
from ._core import Molecule, MoleculeCollection
from ._errors import (
    CalculationCancelledError,
    CalculationError,
    NoExecutablePlanError,
    NumericalFailureError,
)
from ._methods import (
    ExecutionIssue,
    Method,
    MethodCatalog,
    PrerequisiteIssue,
)
from ._parameters import ParameterSet, ParameterSetCatalog
from ._payloads import (
    CalculationProgressPayload,
    ChargeIssuePayload,
    EffectiveCalculationPayload,
    ExecutionIssuePayload,
    ExecutionPolicyPayload,
    FixedIonsPayload,
    InputIssuePayload,
    PrerequisiteIssuePayload,
    RejectionPayload,
)
from ._types import (
    CalculationPhase,
    ChargeIssueKind,
    ExecutionMode,
    ExecutionStatus,
    InputIssueKind,
    MethodOptionValue,
    ParameterMatching,
)
from ._warnings import (
    CHARGE_WARNING_KINDS,
    INPUT_WARNING_KINDS,
    ChargeFWWarning,
    ChargeWarning,
    InputWarning,
    issue_warnings,
)


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class InputIssue:
    """A caution about an input molecule, such as hydrogens that were probably not added."""

    kind: InputIssueKind
    message: str
    molecule_index: int
    conformer_index: int | None = None
    atom_index: int | None = None
    # Every atom the issue concerns, in atom order; atom_index is the one the message names.
    atom_indices: NDArray[np.int64] = field(default_factory=index_array)


def _input_issues(payload: list[InputIssuePayload]) -> tuple[InputIssue, ...]:
    return tuple(
        InputIssue(
            kind=issue["kind"],
            message=issue["message"],
            molecule_index=issue["molecule_index"],
            conformer_index=issue["conformer_index"],
            atom_index=issue["atom_index"],
            atom_indices=index_array(issue["atom_indices"]),
        )
        for issue in payload
    )


@dataclass(frozen=True, slots=True, eq=False, kw_only=True)
class ChargeIssue:
    """A caution about calculated charges: they exist, but are probably wrong."""

    kind: ChargeIssueKind
    message: str
    molecule_index: int
    conformer_index: int | None
    atom_index: int | None
    # Every atom the issue concerns, in atom order; atom_index is the one the message names.
    atom_indices: NDArray[np.int64]


def _charge_issues(payload: list[ChargeIssuePayload]) -> tuple[ChargeIssue, ...]:
    return tuple(
        ChargeIssue(
            kind=issue["kind"],
            message=issue["message"],
            molecule_index=issue["molecule_index"],
            conformer_index=issue["conformer_index"],
            atom_index=issue["atom_index"],
            atom_indices=index_array(issue["atom_indices"]),
        )
        for issue in payload
    )


@dataclass(frozen=True, slots=True, kw_only=True)
class Execution:
    mode: ExecutionMode
    radius: float | None


@dataclass(frozen=True, slots=True, kw_only=True)
class ResourceEstimate:
    """Approximate peak working memory and operation count of executing a plan."""

    memory_bytes: float
    operations: float


@dataclass(frozen=True, slots=True, kw_only=True)
class CalculationProgress:
    """Owned snapshot of one native calculation progress event."""

    phase: CalculationPhase
    mode: ExecutionMode
    method: Method
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
    __slots__ = ("_observer", "_methods")

    def __init__(self, observer: CalculationObserver, methods: MethodCatalog) -> None:
        self._observer = observer
        self._methods = methods

    def on_progress(self, progress: CalculationProgressPayload) -> None:
        self._observer.on_progress(
            CalculationProgress(
                phase=progress["phase"],
                mode=progress["mode"],
                method=self._methods[progress["method_id"]],
                target_index=progress["target_index"],
                target_count=progress["target_count"],
                completed_fragment_count=progress["completed_fragment_count"],
                fragment_count=progress["fragment_count"],
                molecule_index=progress["molecule_index"],
                conformer_index=progress["conformer_index"],
                elapsed_seconds=progress["elapsed_seconds"],
            )
        )

    def cancelled(self) -> bool:
        return self._observer.cancelled()


def _observer_adapter(
    observer: CalculationObserver | None, methods: MethodCatalog
) -> _NativeObserverAdapter | None:
    if observer is None:
        return None
    if not isinstance(observer, CalculationObserver):
        raise TypeError("observer must be a CalculationObserver or None")
    return _NativeObserverAdapter(observer, methods)


def _normalized_threads(threads: int | None) -> int | None:
    return None if threads is None else RequestedCalculation(threads=threads).threads


class CalculationSettings:
    """The choices of a plan, to calculate other molecules with exactly those choices.

    Get it from ``plan.settings`` or ``result.plan.settings`` and pass it as ``settings`` to
    ``calculate()`` or ``assess()``.
    """

    __slots__ = (
        "_method",
        "_parameter_set",
        "_method_options",
        "_parameter_matching",
        "_execution",
        "_fixed_ions",
    )

    _method: Method
    _parameter_set: ParameterSet | None
    _method_options: Mapping[str, MethodOptionValue]
    _parameter_matching: ParameterMatching
    _execution: Execution
    _fixed_ions: tuple[str, ...]

    def __init__(self, *args: object, **kwargs: object) -> None:
        raise TypeError(
            "CalculationSettings objects come from plan.settings and result.plan.settings"
        )

    @classmethod
    def _create(
        cls,
        method: Method,
        parameter_set: ParameterSet | None,
        method_options: Mapping[str, MethodOptionValue],
        execution: Execution,
        requested: RequestedCalculation,
    ) -> CalculationSettings:
        settings = object.__new__(cls)
        settings._method = method
        settings._parameter_set = parameter_set
        settings._method_options = MappingProxyType(dict(method_options))
        settings._parameter_matching = requested.parameter_matching
        settings._execution = execution
        # The requested names, so that ions absent from the assessed molecules are fixed in others.
        settings._fixed_ions = requested.fixed_ions
        return settings

    @property
    def method(self) -> Method:
        return self._method

    @property
    def parameter_set(self) -> ParameterSet | None:
        return self._parameter_set

    @property
    def method_options(self) -> Mapping[str, MethodOptionValue]:
        return self._method_options

    @property
    def parameter_matching(self) -> ParameterMatching:
        return self._parameter_matching

    @property
    def execution(self) -> Execution:
        return self._execution

    @property
    def fixed_ions(self) -> tuple[str, ...]:
        return self._fixed_ions

    def __eq__(self, other: object) -> bool:
        # Not hashable, because a Method is not.
        if not isinstance(other, CalculationSettings):
            return NotImplemented
        return (
            self._method is other._method
            and self._parameter_set == other._parameter_set
            and self._method_options == other._method_options
            and self._parameter_matching == other._parameter_matching
            and self._execution == other._execution
            and self._fixed_ions == other._fixed_ions
        )

    def __repr__(self) -> str:
        parameter_set = self.parameter_set.id if self.parameter_set is not None else None
        return (
            f"{type(self).__name__}(method={self.method.id!r}, "
            f"parameter_set={parameter_set!r}, method_options={dict(self.method_options)!r}, "
            f"parameter_matching={self.parameter_matching!r}, execution={self.execution!r}, "
            f"fixed_ions={self.fixed_ions!r})"
        )


class Plan:
    """A complete reusable execution plan bound to one assessed molecule collection."""

    __slots__ = (
        "_method",
        "_parameter_set",
        "_method_options",
        "_execution",
        "_warnings",
        "_native",
        "_molecules",
        "_requested",
        "_methods",
        "_parameter_sets",
        "_settings",
    )

    _method: Method
    _parameter_set: ParameterSet | None
    _method_options: Mapping[str, MethodOptionValue]
    _execution: Execution
    _warnings: tuple[ExecutionIssue, ...]
    _native: _native_calculation._NativePlan
    _molecules: MoleculeCollection
    _requested: RequestedCalculation
    _methods: MethodCatalog
    _parameter_sets: ParameterSetCatalog
    _settings: CalculationSettings

    def __init__(self, *args: object, **kwargs: object) -> None:
        raise TypeError("Plan objects come from chargefw.assess()")

    @classmethod
    def _create(
        cls,
        native: _native_calculation._NativePlan,
        molecules: MoleculeCollection,
        methods: MethodCatalog,
        parameter_sets: ParameterSetCatalog,
        requested: RequestedCalculation,
    ) -> Plan:
        plan = object.__new__(cls)
        payload = native.report()
        plan._method = methods[payload["method_id"]]
        plan._parameter_set = _resolved_parameter_set(payload["parameter_set_id"], parameter_sets)
        plan._method_options = MappingProxyType(dict(payload["method_options"]))
        plan._execution = _execution(payload["execution_policy"])
        plan._warnings = tuple(_execution_issue(issue) for issue in payload["warnings"])
        plan._native = native
        plan._molecules = molecules
        plan._requested = requested
        plan._methods = methods
        # Resolves the parameter set of results, including a set that is not bundled.
        plan._parameter_sets = parameter_sets
        plan._settings = CalculationSettings._create(
            plan._method, plan._parameter_set, plan._method_options, plan._execution, requested
        )
        return plan

    @property
    def method(self) -> Method:
        return self._method

    @property
    def parameter_set(self) -> ParameterSet | None:
        return self._parameter_set

    @property
    def method_options(self) -> Mapping[str, MethodOptionValue]:
        return self._method_options

    @property
    def execution(self) -> Execution:
        return self._execution

    @property
    def warnings(self) -> tuple[ExecutionIssue, ...]:
        return self._warnings

    @property
    def settings(self) -> CalculationSettings:
        """The choices of this plan, to calculate other molecules with them."""

        return self._settings

    def estimate(self, *, threads: int | None = None) -> ResourceEstimate:
        """Estimate execution with a thread limit; None uses the assessment's threads."""

        memory_bytes, operations = self._native.estimate(_normalized_threads(threads))
        return ResourceEstimate(memory_bytes=memory_bytes, operations=operations)

    def calculate(
        self, *, threads: int | None = None, observer: CalculationObserver | None = None
    ) -> CalculationResult:
        """Calculate the assessed molecules; threads=None uses the assessment's threads.

        The molecules are not prepared or matched to parameters again.
        """

        result = self._execute(threads, _observer_adapter(observer, self._methods))
        # Raised here rather than in a helper, so that the traceback ends at this method.
        error, cautions = _outcome(result)
        if error is not None:
            raise error
        for message, category in cautions:
            warnings.warn(message, category, stacklevel=2)
        return result

    def _execute(
        self, threads: int | None, observer: _NativeObserverAdapter | None
    ) -> CalculationResult:
        normalized_threads = _normalized_threads(threads)
        requested = (
            self._requested
            if normalized_threads is None
            else replace(self._requested, threads=normalized_threads)
        )
        return CalculationResult._create(
            self._native.calculate(normalized_threads, observer),
            self._molecules,
            requested,
            self._methods,
            self._parameter_sets,
        )

    def __repr__(self) -> str:
        parameter_set = self.parameter_set.id if self.parameter_set is not None else None
        return (
            f"{type(self).__name__}(method={self.method.id!r}, "
            f"parameter_set={parameter_set!r}, mode={self.execution.mode!r})"
        )


@dataclass(frozen=True, slots=True, kw_only=True)
class Rejection:
    method: Method
    parameter_set: ParameterSet | None
    execution: Execution | None
    issues: tuple[PrerequisiteIssue | ExecutionIssue, ...]


@dataclass(frozen=True, slots=True, init=False)
class FixedIons:
    sources: tuple[FixedAtomCharge, ...]

    def __init__(
        self,
        sources: Iterable[FixedAtomCharge],
    ) -> None:
        normalized_sources = tuple(sources)
        if any(not isinstance(source, FixedAtomCharge) for source in normalized_sources):
            raise TypeError("sources must contain only FixedAtomCharge values")
        object.__setattr__(self, "sources", normalized_sources)


@dataclass(frozen=True, slots=True, kw_only=True)
class ExecutedPlan:
    method: Method
    parameter_set: ParameterSet | None
    method_options: Mapping[str, MethodOptionValue]
    execution: Execution
    warnings: tuple[ExecutionIssue, ...]
    fixed_ions: FixedIons | None
    # The choices of this plan, to calculate other molecules with them.
    settings: CalculationSettings = field(repr=False, compare=False)

    def __post_init__(self) -> None:
        object.__setattr__(self, "method_options", MappingProxyType(dict(self.method_options)))


@dataclass(frozen=True, slots=True, kw_only=True)
class CalculationTimings:
    applicability_seconds: float
    computation_seconds: float


def _prerequisite_issue(value: PrerequisiteIssuePayload) -> PrerequisiteIssue:
    return PrerequisiteIssue(
        kind=value["kind"],
        message=value["message"],
        molecule_index=value["molecule_index"],
        atom_index=value["atom_index"],
        bond_index=value["bond_index"],
        conformer_index=value["conformer_index"],
        atom_indices=index_array(value["atom_indices"]),
        bond_indices=index_array(value["bond_indices"]),
    )


def _execution_issue(value: ExecutionIssuePayload) -> ExecutionIssue:
    return ExecutionIssue(
        kind=value["kind"],
        message=value["message"],
        molecule_index=value["molecule_index"],
    )


def _resolved_parameter_set(
    parameter_set_id: str | None, parameter_sets: ParameterSetCatalog
) -> ParameterSet | None:
    return None if parameter_set_id is None else parameter_sets[parameter_set_id]


def _rejections(
    values: list[RejectionPayload],
    methods: MethodCatalog,
    parameter_sets: ParameterSetCatalog,
) -> tuple[Rejection, ...]:
    return tuple(
        Rejection(
            method=methods[value["method_id"]],
            parameter_set=_resolved_parameter_set(value["parameter_set_id"], parameter_sets),
            execution=(
                None if value["execution_policy"] is None else _execution(value["execution_policy"])
            ),
            issues=(
                *(_prerequisite_issue(issue) for issue in value["prerequisite_issues"]),
                *(_execution_issue(issue) for issue in value["execution_issues"]),
            ),
        )
        for value in values
    )


def _execution(value: ExecutionPolicyPayload) -> Execution:
    return Execution(
        mode=value["mode"],
        radius=value["radius"],
    )


def _component_id(molecule: Molecule, atom_index: int) -> str | None:
    # The residue that fixed-ion selection matched, so the name is the one the caller requested.
    mapping = molecule.source_mapping
    if mapping is None:
        return None
    for component in mapping.components:
        if atom_index in component.atom_indices:
            return component.component_id
    return None


def _fixed_ions_provenance(
    value: FixedIonsPayload | None, molecules: MoleculeCollection
) -> FixedIons | None:
    if value is None:
        return None
    return FixedIons(
        sources=tuple(
            FixedAtomCharge(
                **source,
                component_id=_component_id(
                    molecules[source["molecule_index"]], source["atom_index"]
                ),
            )
            for source in value["sources"]
        ),
    )


def _executed_plan(
    value: EffectiveCalculationPayload | None,
    methods: MethodCatalog,
    parameter_sets: ParameterSetCatalog,
    molecules: MoleculeCollection,
    requested: RequestedCalculation,
) -> ExecutedPlan | None:
    if value is None:
        return None
    method = methods[value["method_id"]]
    parameter_set = _resolved_parameter_set(value["parameter_set_id"], parameter_sets)
    execution = _execution(value["execution_policy"])
    return ExecutedPlan(
        method=method,
        parameter_set=parameter_set,
        method_options=value["method_options"],
        execution=execution,
        warnings=tuple(_execution_issue(issue) for issue in value["execution_issues"]),
        fixed_ions=_fixed_ions_provenance(value["fixed_ions"], molecules),
        settings=CalculationSettings._create(
            method, parameter_set, value["method_options"], execution, requested
        ),
    )


_REJECTIONS_IN_MESSAGE = 5


def _no_executable_plan_message(rejections: tuple[Rejection, ...]) -> str:
    """List the first reason of the first rejected candidates; the result holds all of them."""

    lines = ["no method and parameter set apply to the molecules:"]
    for rejection in rejections[:_REJECTIONS_IN_MESSAGE]:
        candidate = f"method '{rejection.method.id}'"
        if rejection.parameter_set is not None:
            candidate += f", parameter set '{rejection.parameter_set.id}'"
        line = f"  {candidate}"
        if rejection.issues:
            line += f": {rejection.issues[0].message}"
        if len(rejection.issues) > 1:
            line += f" (and {len(rejection.issues) - 1} more)"
        lines.append(line)
    if len(rejections) > _REJECTIONS_IN_MESSAGE:
        omitted = len(rejections) - _REJECTIONS_IN_MESSAGE
        lines.append(f"  ... {omitted} more rejected candidates in error.result.rejections")
    return "\n".join(lines)


class CalculationResult:
    """Owned calculation status, assignments, diagnostics, and provenance."""

    __slots__ = (
        "_status",
        "_native",
        "_molecules",
        "_requested",
        "_assignments",
        "_assignments_by_molecule",
        "_rejections",
        "_input_warnings",
        "_charge_warnings",
        "_plan",
        "_failure_message",
        "_timings",
    )

    _status: ExecutionStatus
    _native: _native_calculation._NativeExecutionResult
    _molecules: MoleculeCollection
    _requested: RequestedCalculation
    _assignments: tuple[ChargeAssignment, ...]
    _assignments_by_molecule: tuple[tuple[ChargeAssignment, ...], ...]
    _rejections: tuple[Rejection, ...]
    _input_warnings: tuple[InputIssue, ...]
    _charge_warnings: tuple[ChargeIssue, ...]
    _plan: ExecutedPlan | None
    _failure_message: str | None
    _timings: CalculationTimings

    def __init__(self, *args: object, **kwargs: object) -> None:
        raise TypeError(
            "CalculationResult objects come from chargefw.calculate() and Plan.calculate()"
        )

    @classmethod
    def _create(
        cls,
        native: _native_calculation._NativeExecutionResult,
        molecules: MoleculeCollection,
        requested: RequestedCalculation,
        methods: MethodCatalog,
        parameter_sets: ParameterSetCatalog,
        rejections: tuple[Rejection, ...] = (),
    ) -> CalculationResult:
        result = object.__new__(cls)
        payload = native.report()
        assignments: list[ChargeAssignment] = []
        assignments_by_molecule: list[list[ChargeAssignment]] = [[] for _ in molecules]
        charges = payload["charges"]
        if charges is not None:
            for item in charges["assignments"]:
                molecule_index = int(item["molecule_index"])
                conformer_index = item["conformer_index"]
                molecule = molecules[molecule_index]
                assignment = ChargeAssignment._from_native_values(
                    item["values"],
                    molecule_index=molecule_index,
                    conformer_index=(int(conformer_index) if conformer_index is not None else None),
                    source=molecule.source,
                    atom_ids=molecule.atom_ids,
                )
                assignments.append(assignment)
                assignments_by_molecule[molecule_index].append(assignment)
        result._native = native
        result._molecules = molecules
        result._status = payload["status"]
        result._requested = requested
        result._assignments = tuple(assignments)
        result._assignments_by_molecule = tuple(tuple(group) for group in assignments_by_molecule)
        result._rejections = rejections or _rejections(
            payload["rejections"], methods, parameter_sets
        )
        result._input_warnings = _input_issues(payload["input_warnings"])
        result._charge_warnings = _charge_issues(payload["charge_warnings"])
        result._plan = _executed_plan(
            payload["effective"], methods, parameter_sets, molecules, requested
        )
        result._failure_message = payload["failure_message"]
        result._timings = CalculationTimings(
            applicability_seconds=payload["metrics"]["applicability_seconds"],
            computation_seconds=payload["metrics"]["computation_seconds"],
        )
        return result

    @property
    def status(self) -> ExecutionStatus:
        return self._status

    @property
    def molecules(self) -> MoleculeCollection:
        return self._molecules

    @property
    def assignments(self) -> tuple[ChargeAssignment, ...]:
        return self._assignments

    @property
    def assignments_by_molecule(self) -> tuple[tuple[ChargeAssignment, ...], ...]:
        """Return source-ordered assignment groups aligned with the input molecule collection."""

        return self._assignments_by_molecule

    def assignment(
        self, *, molecule_index: int, conformer_index: int | None = None
    ) -> ChargeAssignment:
        """Return charges for one source molecule and, when applicable, one conformer."""

        molecule_index = as_integer(molecule_index, "molecule_index")
        if molecule_index < 0 or molecule_index >= len(self._assignments_by_molecule):
            raise IndexError("molecule_index is outside the calculation input collection")

        assignments = self._assignments_by_molecule[molecule_index]
        if conformer_index is None:
            matches = tuple(item for item in assignments if item.conformer_index is None)
            if len(matches) == 1:
                return matches[0]
            if assignments:
                raise ValueError("conformer_index is required for conformer-specific assignments")
            raise IndexError("molecule has no charge assignment")

        conformer_index = as_integer(
            conformer_index, "conformer_index", expected="an integer or None"
        )
        if conformer_index < 0:
            raise IndexError("conformer_index must be non-negative")
        matches = tuple(item for item in assignments if item.conformer_index == conformer_index)
        if len(matches) == 1:
            return matches[0]
        if any(item.conformer_index is None for item in assignments):
            raise ValueError("molecule has no conformer-specific assignments")
        raise IndexError("conformer has no charge assignment")

    @property
    def rejections(self) -> tuple[Rejection, ...]:
        return self._rejections

    @property
    def requested(self) -> RequestedCalculation:
        return self._requested

    @property
    def plan(self) -> ExecutedPlan | None:
        return self._plan

    @property
    def input_warnings(self) -> tuple[InputIssue, ...]:
        """Cautions about the input molecules, such as hydrogens that were probably not added."""

        return self._input_warnings

    @property
    def charge_warnings(self) -> tuple[ChargeIssue, ...]:
        """Cautions about the calculated charges, such as charges far from the formal charges."""

        return self._charge_warnings

    @property
    def failure_message(self) -> str | None:
        return self._failure_message

    @property
    def timings(self) -> CalculationTimings:
        return self._timings

    def _status_error(self) -> CalculationError | None:
        exception_types = {
            "no_executable_plan": NoExecutablePlanError,
            "numerical_failure": NumericalFailureError,
            "cancelled": CalculationCancelledError,
        }
        exception_type = exception_types.get(self.status)
        if exception_type is None:
            return None
        if self.status == "no_executable_plan":
            return exception_type(_no_executable_plan_message(self.rejections), self)
        return exception_type(self.failure_message or self.status, self)


def _outcome(
    result: CalculationResult,
) -> tuple[CalculationError | None, list[tuple[str, type[ChargeFWWarning]]]]:
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
