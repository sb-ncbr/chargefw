"""Exceptions that ChargeFW raises for input it cannot read or calculations that fail.

Wrong arguments raise TypeError or ValueError, and files that cannot be opened OSError, as elsewhere
in Python.
"""

from __future__ import annotations

from collections.abc import Iterable
from functools import partial
from typing import TYPE_CHECKING, Any

from ._arrays import index_array
from ._public import set_public_module
from ._types import OutputErrorKind, ParseErrorKind

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
    kinds may be added. ``reason`` is why as text; the message of the error says where before it.
    ``source_name`` names the source, or is empty.
    Where the reader knows them, ``line`` is the line of the source, counted from one,
    ``record_index`` the record, counted from zero, and ``record_id`` the record's identifier in the
    source, such as its title; otherwise they are ``None``. ``line_text`` is the text of that line,
    without its line ending, or ``None`` when there is no line, the input ends before it, or the
    file could not be decoded.
    """

    def __init__(
        self,
        reason: str,
        *,
        kind: ParseErrorKind,
        source_name: str = "",
        line: int | None = None,
        line_text: str | None = None,
        record_index: int | None = None,
        record_id: str | int | None = None,
    ) -> None:
        # The message says where, then why: "input.sdf:24: record 3 ('title'): reason".
        where = ""
        if line is not None:
            where = f"{source_name}:{line}: " if source_name else f"line {line}: "
        elif source_name:
            where = f"{source_name}: "
        if record_index is not None:
            name = "" if record_id in (None, "") else f" ('{record_id}')"
            where += f"record {record_index + 1}{name}: "
        super().__init__(where + reason)
        self.kind = kind
        self.reason = reason
        self.source_name = source_name
        self.line = line
        self.line_text = line_text
        self.record_index = record_index
        self.record_id = record_id

    def __reduce__(self) -> tuple[Any, ...]:
        # Rebuilt from the reason and the fields, which give the message; the state restores the
        # rest, such as notes and attributes the caller added.
        rebuild = partial(
            type(self),
            kind=self.kind,
            source_name=self.source_name,
            line=self.line,
            line_text=self.line_text,
            record_index=self.record_index,
            record_id=self.record_id,
        )
        return rebuild, (self.reason,), vars(self)


class OutputError(ChargeFWError, ValueError):
    """A result cannot be written in the requested format, or charges cannot be added to a target.

    ``kind`` is the cause as a code: ``"missing_coordinates"`` for a molecule without coordinates in
    a format that needs them, ``"unsupported_shape"`` for a number of molecules, conformers, or
    atoms the format cannot hold, ``"unrepresentable"`` for a value the format cannot hold,
    ``"target_mismatch"`` for a document or toolkit molecule that is not the one that was read, and
    ``"existing_charges"`` for charges already in the target that block the request. More kinds may
    be added. ``reason`` is why as text; the message of the error says where before it.

    ``molecule_index`` is the input molecule the error concerns, or ``None`` when it concerns the
    result as a whole. The other fields say where in the molecule, as for an issue:
    ``conformer_index``, ``atom_index``, ``other_atom_index`` for the second atom of a pair,
    ``bond_index``, and the arrays ``atom_indices`` and ``bond_indices``. All count from zero.
    """

    def __init__(
        self,
        reason: str,
        *,
        kind: OutputErrorKind,
        molecule_index: int | None = None,
        conformer_index: int | None = None,
        atom_index: int | None = None,
        other_atom_index: int | None = None,
        bond_index: int | None = None,
        atom_indices: Iterable[int] = (),
        bond_indices: Iterable[int] = (),
    ) -> None:
        # The message says where, counted from one, then why.
        where = [
            f"{name} {index + 1}"
            for name, index in (("molecule", molecule_index), ("conformer", conformer_index))
            if index is not None
        ]
        if atom_index is not None and other_atom_index is not None:
            where.append(f"atoms {atom_index + 1} and {other_atom_index + 1}")
        elif atom_index is not None:
            where.append(f"atom {atom_index + 1}")
        if bond_index is not None:
            where.append(f"bond {bond_index + 1}")
        super().__init__(f"{', '.join(where)}: {reason}" if where else reason)
        self.kind = kind
        self.reason = reason
        self.molecule_index = molecule_index
        self.conformer_index = conformer_index
        self.atom_index = atom_index
        self.other_atom_index = other_atom_index
        self.bond_index = bond_index
        self.atom_indices = index_array(atom_indices)
        self.bond_indices = index_array(bond_indices)

    def __reduce__(self) -> tuple[Any, ...]:
        # Rebuilt from the reason and the fields, which give the message; the state restores the
        # rest, such as notes.
        rebuild = partial(
            type(self),
            kind=self.kind,
            molecule_index=self.molecule_index,
            conformer_index=self.conformer_index,
            atom_index=self.atom_index,
            other_atom_index=self.other_atom_index,
            bond_index=self.bond_index,
        )
        return rebuild, (self.reason,), vars(self)

    def __setstate__(self, state: dict[str, Any] | None) -> None:
        vars(self).update(state or {})
        # Unpickled arrays are writable; these are not.
        self.atom_indices = index_array(self.atom_indices)
        self.bond_indices = index_array(self.bond_indices)


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
    OutputError,
    CalculationError,
    NoExecutablePlanError,
    NumericalFailureError,
    CalculationCancelledError,
)
