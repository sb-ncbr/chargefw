"""Fitting parameter values of ChargeFW methods to reference charges."""

from ._diagnostics import ChargeDiagnostics, Diagnostics
from ._evaluator import EvaluationError, Evaluator, Slot
from ._fit import FitResult, fit
from ._observer import FitObserver, FitProgress
from ._template import make_parameter_template

for _value_type in (
    ChargeDiagnostics,
    Diagnostics,
    EvaluationError,
    Evaluator,
    FitObserver,
    FitProgress,
    FitResult,
    Slot,
):
    _value_type.__module__ = __name__

__all__ = [
    "ChargeDiagnostics",
    "Diagnostics",
    "EvaluationError",
    "Evaluator",
    "FitObserver",
    "FitProgress",
    "FitResult",
    "Slot",
    "fit",
    "make_parameter_template",
]
