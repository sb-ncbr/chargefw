"""Shared private type aliases."""

from typing import Literal

type MethodOptionValue = bool | int | float | str
type PortableId = str | int
type RequestedExecution = Literal["auto", "full", "cutoff", "cover"]
type ParameterMatching = Literal["strict", "permissive"]
type MethodScope = Literal["recommended", "all"]
type AtomClassifier = Literal["plain", "hbo", "bonded"]
type BondClassifier = Literal["plain", "bo"]
type ExecutionMode = Literal["full", "cutoff", "cover"]
type CalculationPhase = Literal[
    "computation_started",
    "computation_finished",
    "target_started",
    "target_finished",
    "fragment_progress",
]
type ExecutionIssueKind = Literal["resource_threshold_exceeded", "unsupported_execution_mode"]
type InputIssueKind = Literal["hydrogens_missing", "atoms_too_close"]
type ImportIssueKind = Literal[
    "partial_charges_ignored",
    "coordinates_2d_ignored",
    "v2000_property_ignored",
    "v3000_bond_order_ignored",
    "aromatic_bonds_converted",
    "symmetry_connections_ignored",
    "unbonded_residue_atoms",
]
type PrerequisiteIssueKind = Literal[
    "invalid_options",
    "missing_feature",
    "invalid_geometry",
    "unsupported_molecule",
    "missing_parameters",
    "parameter_classification_failed",
    "unsupported_fixed_ions",
]
type MethodOptionType = Literal["boolean", "integer", "floating_point", "string"]
type ExecutionStatus = Literal[
    "success",
    "no_executable_plan",
    "numerical_failure",
    "cancelled",
]
