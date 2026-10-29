"""ChargeFW Python API."""

from typing import Final

from . import io as io
from ._chargefw import version as _native_version
from ._chargefw.calculation import _fixed_ion_names
from ._elements import Element, PeriodicTable, periodic_table
from ._methods import ExecutionIssue, Method, MethodOption, PrerequisiteIssue
from ._parameters import AtomParameterKey, BondParameterKey, ParameterSet
from .calculation import (
    Assessment,
    CalculationCancelledError,
    CalculationObserver,
    CalculationProgress,
    CalculationResult,
    CalculationTimings,
    ChargeFWError,
    ExecutedPlan,
    ExecutionPolicy,
    FixedAtomCharge,
    FixedIons,
    NoExecutablePlanError,
    NumericalFailureError,
    Plan,
    Rejection,
    RequestedCalculation,
    ResourceEstimate,
    assess,
    calculate,
    methods,
    parameter_sets,
)
from .charges import ChargeAssignment
from .core import (
    ImportDiagnostic,
    Molecule,
    MoleculeCollection,
    PortableId,
    SourceAtomReference,
    SourceComponentInstance,
    SourceConformerReference,
    SourceHierarchyLabels,
    SourceIdentity,
    SourceMapping,
    SourceStructuralLabels,
)

__version__ = _native_version()
COMMON_IONS: Final[tuple[str, ...]] = _fixed_ion_names(True)
ALL_IONS: Final[tuple[str, ...]] = _fixed_ion_names(False)

__all__ = [
    "__version__",
    "COMMON_IONS",
    "ALL_IONS",
    "Molecule",
    "MoleculeCollection",
    "PortableId",
    "SourceIdentity",
    "SourceAtomReference",
    "SourceComponentInstance",
    "SourceConformerReference",
    "SourceHierarchyLabels",
    "SourceStructuralLabels",
    "SourceMapping",
    "ImportDiagnostic",
    "ChargeAssignment",
    "FixedAtomCharge",
    "CalculationResult",
    "Assessment",
    "CalculationObserver",
    "CalculationProgress",
    "assess",
    "calculate",
    "methods",
    "parameter_sets",
    "periodic_table",
    "RequestedCalculation",
    "ChargeFWError",
    "NoExecutablePlanError",
    "NumericalFailureError",
    "CalculationCancelledError",
    "PrerequisiteIssue",
    "ExecutionIssue",
    "ExecutionPolicy",
    "Plan",
    "Rejection",
    "ResourceEstimate",
    "ExecutedPlan",
    "FixedIons",
    "CalculationTimings",
    "MethodOption",
    "Method",
    "ParameterSet",
    "Element",
    "PeriodicTable",
    "AtomParameterKey",
    "BondParameterKey",
    "io",
]
