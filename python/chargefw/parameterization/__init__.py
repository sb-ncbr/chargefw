"""Fitting parameter values of ChargeFW methods to reference charges."""

from .._public import set_public_module as _set_public_module
from ._diagnostics import ChargeDiagnostics, Diagnostics
from ._evaluator import EvaluationError, Evaluator, Slot
from ._fit import FitResult, fit
from ._observer import FitObserver, FitProgress
from ._template import make_parameter_template

_set_public_module(
    __name__,
    ChargeDiagnostics,
    Diagnostics,
    EvaluationError,
    Evaluator,
    FitObserver,
    FitProgress,
    FitResult,
    Slot,
)

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
