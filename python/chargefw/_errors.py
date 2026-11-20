"""Exceptions that ChargeFW raises for input it cannot read or calculations that fail.

Wrong arguments raise TypeError or ValueError, and files that cannot be opened OSError, as elsewhere
in Python.
"""

from __future__ import annotations

from functools import partial
from typing import TYPE_CHECKING, Any

from ._public import set_public_module
from ._types import ParseErrorKind

if TYPE_CHECKING:
    from ._calculation_values import CalculationResult


class ChargeFWError(Exception):
    """Base class of the errors ChargeFW raises."""


class ParseError(ChargeFWError, ValueError):
    """Molecular input or a parameter-set file cannot be read.

    ``kind`` is the cause as a code: ``"no_records"`` for input with nothing to read, such as an
    empty file or a file of another format, ``"invalid_encoding"`` for text that cannot be decoded,
    ``"malformed"`` for input that breaks its format, ``"unsupported"`` for valid data that ChargeFW
    does not read, and ``"inconsistent_models"`` for PDB or mmCIF models with different atoms. More
    kinds may be added. ``reason`` is why as text, without where ChargeFW knows. ``source_name`` names the source, or is empty.
    Where the reader knows them, ``line`` is the line of the source, counted from one,
    ``record_index`` the record, counted from zero, and ``record_id`` the record's identifier in the
    source, such as its title; otherwise they are ``None``. ``line_text`` is the text of that line,
    without its line ending, or ``None`` when there is no line, the input ends before it, or the
    file could not be decoded.
    """

    def __init__(
        self,
        message: str,
        *,
        kind: ParseErrorKind,
        reason: str | None = None,
        source_name: str = "",
        line: int | None = None,
        line_text: str | None = None,
        record_index: int | None = None,
        record_id: str | int | None = None,
    ) -> None:
        super().__init__(message)
        self.kind = kind
        self.reason = message if reason is None else reason
        self.source_name = source_name
        self.line = line
        self.line_text = line_text
        self.record_index = record_index
        self.record_id = record_id

    def __reduce__(self) -> tuple[Any, ...]:
        # Exceptions are rebuilt from their positional arguments alone, which lacks the required
        # kind. The state restores the rest, such as notes and attributes the caller added.
        rebuild = partial(
            type(self),
            kind=self.kind,
            reason=self.reason,
            source_name=self.source_name,
            line=self.line,
            line_text=self.line_text,
            record_index=self.record_index,
            record_id=self.record_id,
        )
        return rebuild, self.args, vars(self)


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
