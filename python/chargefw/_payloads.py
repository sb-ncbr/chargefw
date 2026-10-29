"""Typed private payloads returned by the nanobind extension."""

from __future__ import annotations

from typing import TypedDict

import numpy as np
from numpy.typing import NDArray

from ._types import (
    AtomClassifier,
    BondClassifier,
    CalculationPhase,
    ExecutionIssueKind,
    ExecutionMode,
    ExecutionStatus,
    MethodOptionType,
    MethodOptionValue,
    PrerequisiteIssueKind,
)


class CalculationProgressPayload(TypedDict):
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


class PrerequisiteIssuePayload(TypedDict):
    kind: PrerequisiteIssueKind
    message: str
    molecule_index: int | None
    atom_index: int | None
    bond_index: int | None
    conformer_index: int | None


class ExecutionIssuePayload(TypedDict):
    kind: ExecutionIssueKind
    message: str
    molecule_index: int | None


class ExecutionPolicyPayload(TypedDict):
    mode: ExecutionMode
    radius: float | None


class FixedAtomChargePayload(TypedDict):
    molecule_index: int
    atom_index: int
    charge: float


class FixedIonsPayload(TypedDict):
    sources: list[FixedAtomChargePayload]


class EffectiveCalculationPayload(TypedDict):
    method_id: str
    parameter_set_id: str | None
    method_options: dict[str, MethodOptionValue]
    execution_policy: ExecutionPolicyPayload
    execution_issues: list[ExecutionIssuePayload]
    fixed_ions: FixedIonsPayload | None


class ExecutionPlanPayload(TypedDict):
    method_id: str
    parameter_set_id: str | None
    method_options: dict[str, MethodOptionValue]
    execution_policy: ExecutionPolicyPayload
    warnings: list[ExecutionIssuePayload]


class RejectionPayload(TypedDict):
    method_id: str
    parameter_set_id: str | None
    execution_policy: ExecutionPolicyPayload | None
    prerequisite_issues: list[PrerequisiteIssuePayload]
    execution_issues: list[ExecutionIssuePayload]


class ChargeAssignmentPayload(TypedDict):
    molecule_index: int
    conformer_index: int | None
    values: NDArray[np.float64]


class ChargeSetPayload(TypedDict):
    assignments: list[ChargeAssignmentPayload]


class CalculationMetricsPayload(TypedDict):
    applicability_seconds: float
    computation_seconds: float


class ExecutionResultPayload(TypedDict):
    status: ExecutionStatus
    rejections: list[RejectionPayload]
    failure_message: str | None
    metrics: CalculationMetricsPayload
    effective: EffectiveCalculationPayload | None
    charges: ChargeSetPayload | None


class AssessmentReportPayload(TypedDict):
    rejections: list[RejectionPayload]
    applicability_seconds: float


class MethodOptionDescriptorPayload(TypedDict):
    id: str
    description: str
    type: MethodOptionType
    default: MethodOptionValue
    choices: list[MethodOptionValue]
    minimum: MethodOptionValue | None
    minimum_inclusive: bool
    maximum: MethodOptionValue | None
    maximum_inclusive: bool


class MethodDescriptorPayload(TypedDict):
    id: str
    name: str
    full_name: str
    publication: str | None
    notes: str
    requires_coordinates: bool
    time_complexity: str
    memory_complexity: str
    supports_cutoff: bool
    supports_cover: bool
    supports_fixed_ions: bool
    options: list[MethodOptionDescriptorPayload]


class ElementPayload(TypedDict):
    atomic_number: int
    symbol: str
    name: str
    period: int
    group: int | None


class AtomParameterKeyPayload(TypedDict):
    element: str
    classifier: AtomClassifier
    type: str


class BondParameterKeyPayload(TypedDict):
    atoms: tuple[AtomParameterKeyPayload, AtomParameterKeyPayload]
    classifier: BondClassifier
    type: str


class ParameterSetDescriptorPayload(TypedDict):
    id: str
    method_id: str
    name: str
    publication: str
    notes: str
    atom_keys: list[AtomParameterKeyPayload]
    bond_keys: list[BondParameterKeyPayload]
