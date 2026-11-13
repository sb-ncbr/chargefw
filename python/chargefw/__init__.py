"""ChargeFW Python API."""

from typing import Final

from . import _cpu as _cpu  # Checks the CPU before the extension is loaded.
from . import io as io
from ._chargefw import version as _native_version
from ._chargefw.calculation import _fixed_ion_names
from ._elements import Element, PeriodicTable, periodic_table
from ._methods import ExecutionIssue, Method, MethodOption, PrerequisiteIssue
from ._parameters import AtomParameterKey, BondParameterKey, ParameterSet
from ._warnings import ChargeFWWarning, ChargeWarning, InputWarning
from .calculation import (
    Assessment,
    CalculationCancelledError,
    CalculationObserver,
    CalculationProgress,
    CalculationResult,
    CalculationSettings,
    CalculationTimings,
    ChargeFWError,
    ChargeIssue,
    ExecutedPlan,
    Execution,
    FixedAtomCharge,
    FixedIons,
    InputIssue,
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
    read_parameter_set,
    write_parameter_set,
)
from .charges import ChargeAssignment
from .core import (
    ImportIssue,
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
    "ImportIssue",
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
    "read_parameter_set",
    "write_parameter_set",
    "RequestedCalculation",
    "ChargeFWError",
    "NoExecutablePlanError",
    "NumericalFailureError",
    "CalculationCancelledError",
    "PrerequisiteIssue",
    "ExecutionIssue",
    "InputIssue",
    "ChargeIssue",
    "ChargeFWWarning",
    "InputWarning",
    "ChargeWarning",
    "Execution",
    "Plan",
    "Rejection",
    "ResourceEstimate",
    "ExecutedPlan",
    "FixedIons",
    "CalculationTimings",
    "CalculationSettings",
    "MethodOption",
    "Method",
    "ParameterSet",
    "Element",
    "PeriodicTable",
    "AtomParameterKey",
    "BondParameterKey",
    "io",
]
