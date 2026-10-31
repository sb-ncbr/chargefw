"""Private immutable calculation result values and native-payload conversion."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from operator import index as as_index
from types import MappingProxyType

from ._calculation_options import FixedAtomCharge, RequestedCalculation
from ._chargefw import calculation as _native_calculation
from ._methods import (
    ExecutionIssue,
    Method,
    MethodCatalog,
    PrerequisiteIssue,
)
from ._parameters import ParameterSet, ParameterSetCatalog
from ._payloads import (
    EffectiveCalculationPayload,
    ExecutionIssuePayload,
    ExecutionPolicyPayload,
    FixedIonsPayload,
    PrerequisiteIssuePayload,
    RejectionPayload,
)
from ._types import ExecutionMode, ExecutionStatus, MethodOptionValue
from .charges import ChargeAssignment
from .core import Molecule, MoleculeCollection


@dataclass(frozen=True, slots=True)
class Execution:
    mode: ExecutionMode
    radius: float | None


@dataclass(frozen=True, slots=True)
class ResourceEstimate:
    """Approximate peak working memory and operation count of executing a plan."""

    memory_bytes: float
    operations: float


class Plan:
    """A complete reusable execution plan bound to one assessed molecule collection."""

    __slots__ = (
        "_method",
        "_parameter_set",
        "_options",
        "_execution",
        "_warnings",
        "_native",
        "_molecules",
        "_requested",
    )

    def __init__(
        self,
        native: _native_calculation._NativePlan,
        molecules: MoleculeCollection,
        methods: MethodCatalog,
        parameter_sets: ParameterSetCatalog,
        requested: RequestedCalculation,
    ) -> None:
        payload = native.report()
        self._method = methods[payload["method_id"]]
        self._parameter_set = _resolved_parameter_set(payload["parameter_set_id"], parameter_sets)
        self._options = MappingProxyType(dict(payload["method_options"]))
        self._execution = _execution(payload["execution_policy"])
        self._warnings = tuple(_execution_issue(issue) for issue in payload["warnings"])
        self._native = native
        self._molecules = molecules
        self._requested = requested

    @property
    def method(self) -> Method:
        return self._method

    @property
    def parameter_set(self) -> ParameterSet | None:
        return self._parameter_set

    @property
    def options(self) -> Mapping[str, MethodOptionValue]:
        return self._options

    @property
    def execution(self) -> Execution:
        return self._execution

    @property
    def warnings(self) -> tuple[ExecutionIssue, ...]:
        return self._warnings

    def _matches(self, molecules: MoleculeCollection) -> bool:
        return (
            molecules.name == self._molecules.name
            and len(molecules) == len(self._molecules)
            and all(
                supplied is assessed
                for supplied, assessed in zip(molecules, self._molecules, strict=True)
            )
        )

    def estimate(self, threads: int | None = None) -> ResourceEstimate:
        """Estimate execution with a thread limit; None uses the assessment's threads."""

        normalized_threads = (
            None if threads is None else RequestedCalculation(threads=threads).threads
        )
        memory_bytes, operations = self._native.estimate(normalized_threads)
        return ResourceEstimate(memory_bytes=memory_bytes, operations=operations)

    def _calculate(
        self, threads: int | None, observer: object | None
    ) -> _native_calculation._NativeExecutionResult:
        return self._native.calculate(threads, observer)

    def __repr__(self) -> str:
        parameter_set = self.parameter_set.id if self.parameter_set is not None else None
        return (
            f"{type(self).__name__}(method={self.method.id!r}, "
            f"parameter_set={parameter_set!r}, mode={self.execution.mode!r})"
        )


@dataclass(frozen=True, slots=True)
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


@dataclass(frozen=True, slots=True)
class ExecutedPlan:
    method: Method
    parameter_set: ParameterSet | None
    options: Mapping[str, MethodOptionValue]
    execution: Execution
    warnings: tuple[ExecutionIssue, ...]
    fixed_ions: FixedIons | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "options", MappingProxyType(dict(self.options)))


@dataclass(frozen=True, slots=True)
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
) -> ExecutedPlan | None:
    if value is None:
        return None
    return ExecutedPlan(
        method=methods[value["method_id"]],
        parameter_set=_resolved_parameter_set(value["parameter_set_id"], parameter_sets),
        options=value["method_options"],
        execution=_execution(value["execution_policy"]),
        warnings=tuple(_execution_issue(issue) for issue in value["execution_issues"]),
        fixed_ions=_fixed_ions_provenance(value["fixed_ions"], molecules),
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


class ChargeFWError(RuntimeError):
    """Base class for a failed ChargeFW calculation result."""

    def __init__(self, message: str, result: CalculationResult) -> None:
        super().__init__(message)
        self.result = result


class NoExecutablePlanError(ChargeFWError):
    pass


class NumericalFailureError(ChargeFWError):
    pass


class CalculationCancelledError(ChargeFWError):
    pass


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
    _plan: ExecutedPlan | None
    _failure_message: str | None
    _timings: CalculationTimings

    def __init__(
        self,
        native: _native_calculation._NativeExecutionResult,
        molecules: MoleculeCollection,
        requested: RequestedCalculation,
        methods: MethodCatalog,
        parameter_sets: ParameterSetCatalog,
        rejections: tuple[Rejection, ...] = (),
    ) -> None:
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
        self._native = native
        self._molecules = molecules
        self._status = payload["status"]
        self._requested = requested
        self._assignments = tuple(assignments)
        self._assignments_by_molecule = tuple(tuple(group) for group in assignments_by_molecule)
        self._rejections = rejections or _rejections(payload["rejections"], methods, parameter_sets)
        self._plan = _executed_plan(payload["effective"], methods, parameter_sets, molecules)
        self._failure_message = payload["failure_message"]
        self._timings = CalculationTimings(
            applicability_seconds=payload["metrics"]["applicability_seconds"],
            computation_seconds=payload["metrics"]["computation_seconds"],
        )

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

    def assignment(self, *, molecule: int, conformer: int | None = None) -> ChargeAssignment:
        """Return charges for one source molecule and, when applicable, one conformer."""

        if isinstance(molecule, bool):
            raise TypeError("molecule must be an integer")
        try:
            molecule_index = as_index(molecule)
        except TypeError as error:
            raise TypeError("molecule must be an integer") from error
        if molecule_index < 0 or molecule_index >= len(self._assignments_by_molecule):
            raise IndexError("molecule index is outside the calculation input collection")

        assignments = self._assignments_by_molecule[molecule_index]
        if conformer is None:
            matches = tuple(item for item in assignments if item.conformer_index is None)
            if len(matches) == 1:
                return matches[0]
            if assignments:
                raise ValueError("conformer is required for conformer-specific assignments")
            raise IndexError("molecule has no charge assignment")

        if isinstance(conformer, bool):
            raise TypeError("conformer must be an integer or None")
        try:
            conformer_index = as_index(conformer)
        except TypeError as error:
            raise TypeError("conformer must be an integer or None") from error
        if conformer_index < 0:
            raise IndexError("conformer index must be non-negative")
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
    def warnings(self) -> tuple[ExecutionIssue, ...]:
        if self._plan is None:
            return ()
        return self._plan.warnings

    @property
    def failure_message(self) -> str | None:
        return self._failure_message

    @property
    def timings(self) -> CalculationTimings:
        return self._timings

    def _raise_for_status(self) -> None:
        exception_types = {
            "no_executable_plan": NoExecutablePlanError,
            "numerical_failure": NumericalFailureError,
            "cancelled": CalculationCancelledError,
        }
        exception_type = exception_types.get(self.status)
        if exception_type is None:
            return
        if self.status == "no_executable_plan":
            raise exception_type(_no_executable_plan_message(self.rejections), self)
        raise exception_type(self.failure_message or self.status, self)
