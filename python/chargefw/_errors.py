"""Exceptions that ChargeFW raises for input it cannot read or calculations that fail.

Wrong arguments raise TypeError or ValueError, and files that cannot be opened OSError, as elsewhere
in Python.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

from ._public import set_public_module

if TYPE_CHECKING:
    from ._calculation_values import CalculationResult


class ChargeFWError(Exception):
    """Base class of the errors ChargeFW raises."""


class ParseError(ChargeFWError, ValueError):
    """Molecular input or a parameter-set file cannot be read."""


class CalculationError(ChargeFWError):
    """A calculation produced no charges; ``result`` holds the failed result."""

    def __init__(self, message: str, result: CalculationResult) -> None:
        super().__init__(message)
        self.result = result


class NoExecutablePlanError(CalculationError):
    """No method and parameter set apply to the molecules with the requested settings."""


class NumericalFailureError(CalculationError):
    """The calculation failed, for example did not converge."""


class CalculationCancelledError(CalculationError):
    """The observer cancelled the calculation."""


set_public_module(
    "chargefw",
    ChargeFWError,
    ParseError,
    CalculationError,
    NoExecutablePlanError,
    NumericalFailureError,
    CalculationCancelledError,
)
