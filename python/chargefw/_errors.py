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
    """Molecular input or a parameter-set file cannot be read.

    ``reason`` is why, without where ChargeFW knows. ``source_name`` names the source, or is empty.
    Where the reader knows them, ``line`` is the line of the source, counted from one,
    ``record_index`` the record, counted from zero, and ``record_id`` the record's identifier in the
    source, such as its title; otherwise they are ``None``.
    """

    def __init__(
        self,
        message: str,
        *,
        reason: str | None = None,
        source_name: str = "",
        line: int | None = None,
        record_index: int | None = None,
        record_id: str | int | None = None,
    ) -> None:
        super().__init__(message)
        self.reason = message if reason is None else reason
        self.source_name = source_name
        self.line = line
        self.record_index = record_index
        self.record_id = record_id


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
