"""Private normalization for Python-facing calculation arguments."""

from __future__ import annotations

import numbers
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from operator import index as as_index
from types import MappingProxyType
from typing import SupportsIndex

import numpy as np

from ._methods import Method
from ._parameters import ParameterSet
from ._types import MethodOptionValue, ParameterMatching, RequestedExecution

_EXECUTIONS = frozenset(("auto", "full", "cutoff", "cover"))


def _normalized_nonnegative_integer(value: SupportsIndex, name: str) -> int:
    if isinstance(value, (bool, np.bool)):
        raise TypeError(f"{name} must be an integer")
    try:
        result = as_index(value)
    except TypeError as error:
        raise TypeError(f"{name} must be an integer") from error
    if result < 0:
        raise ValueError(f"{name} must be non-negative")
    if result > np.iinfo(np.uintp).max:
        raise ValueError(f"{name} is outside the native size range")
    return result


@dataclass(frozen=True, slots=True, init=False)
class FixedAtomCharge:
    """A prescribed fixed charge on one source atom in the input collection."""

    molecule_index: int
    atom_index: int
    charge: float

    def __init__(
        self, molecule_index: SupportsIndex, atom_index: SupportsIndex, charge: float
    ) -> None:
        if isinstance(charge, (bool, np.bool)) or not isinstance(charge, numbers.Real):
            raise TypeError("charge must be a real number")
        try:
            normalized_charge = float(charge)
        except OverflowError as error:
            raise ValueError("charge cannot be represented as a float") from error
        object.__setattr__(
            self,
            "molecule_index",
            _normalized_nonnegative_integer(molecule_index, "molecule_index"),
        )
        object.__setattr__(
            self, "atom_index", _normalized_nonnegative_integer(atom_index, "atom_index")
        )
        object.__setattr__(self, "charge", normalized_charge)


def _normalized_option_values(
    values: Mapping[str, MethodOptionValue], method_id: str
) -> dict[str, MethodOptionValue]:
    if not isinstance(values, Mapping):
        raise TypeError(f"options for method {method_id!r} must be a mapping")
    result: dict[str, MethodOptionValue] = {}
    int_limits = np.iinfo(np.int32)
    for option_id, value in values.items():
        if not isinstance(option_id, str):
            raise TypeError("method option IDs must be strings")
        if isinstance(value, (bool, np.bool)):
            result[option_id] = bool(value)
        elif isinstance(value, numbers.Integral):
            integer = int(value)
            if integer < int_limits.min or integer > int_limits.max:
                raise ValueError(f"method option {method_id}.{option_id} is outside int range")
            result[option_id] = integer
        elif isinstance(value, numbers.Real):
            floating = float(value)
            if not np.isfinite(floating):
                raise ValueError(f"method option {method_id}.{option_id} must be finite")
            result[option_id] = floating
        elif isinstance(value, str):
            result[option_id] = value
        else:
            raise TypeError(
                f"method option {method_id}.{option_id} must be bool, int, float, or str"
            )
    return result


def _normalized_options_by_method(
    values: Mapping[str, Mapping[str, MethodOptionValue]] | None,
) -> dict[str, dict[str, MethodOptionValue]]:
    if values is None:
        return {}
    if not isinstance(values, Mapping):
        raise TypeError("options_by_method must be a mapping")
    result: dict[str, dict[str, MethodOptionValue]] = {}
    for method_id, options in values.items():
        if not isinstance(method_id, str):
            raise TypeError("options_by_method keys must be method ID strings")
        result[method_id] = _normalized_option_values(options, method_id)
    return result


def _frozen_options_by_method(
    values: Mapping[str, Mapping[str, MethodOptionValue]],
) -> Mapping[str, Mapping[str, MethodOptionValue]]:
    return MappingProxyType(
        {method_id: MappingProxyType(dict(options)) for method_id, options in values.items()}
    )


@dataclass(frozen=True, slots=True, init=False)
class RequestedCalculation:
    """Normalized policy requested for an assessment or calculation."""

    method: str | None
    parameter_set: str | None
    options_by_method: Mapping[str, Mapping[str, MethodOptionValue]]
    parameter_matching: ParameterMatching
    execution: RequestedExecution
    radius: float | None
    cutoff_threshold: int | None
    cover_threshold: int | None
    threads: int
    fixed_ions: tuple[str, ...]

    def __init__(
        self,
        *,
        method: str | Method | None = None,
        parameter_set: str | ParameterSet | None = None,
        options: Mapping[str, MethodOptionValue] | None = None,
        options_by_method: Mapping[str, Mapping[str, MethodOptionValue]] | None = None,
        parameter_matching: ParameterMatching = "strict",
        execution: RequestedExecution = "auto",
        radius: float | None = None,
        cutoff_threshold: int | None = 20_000,
        cover_threshold: int | None = 80_000,
        threads: int = 0,
        fixed_ions: Sequence[str] | None = None,
    ) -> None:
        if method is not None and not isinstance(method, (str, Method)):
            raise TypeError("method must be a method ID, Method, or None")
        if parameter_set is not None and not isinstance(parameter_set, (str, ParameterSet)):
            raise TypeError("parameter_set must be a parameter-set ID, ParameterSet, or None")
        method_id = method.id if isinstance(method, Method) else method
        parameter_set_id = (
            parameter_set.id if isinstance(parameter_set, ParameterSet) else parameter_set
        )
        if options is not None and options_by_method is not None:
            raise ValueError("options and options_by_method cannot be used together")
        if options is not None:
            if method_id is None:
                raise ValueError("options requires an explicit method")
            normalized_options = {method_id: _normalized_option_values(options, method_id)}
        else:
            normalized_options = _normalized_options_by_method(options_by_method)

        if not isinstance(parameter_matching, str):
            raise TypeError("parameter_matching must be a string")
        if parameter_matching not in ("strict", "permissive"):
            raise ValueError("parameter_matching must be 'strict' or 'permissive'")
        if not isinstance(execution, str):
            raise TypeError("execution must be a string")
        if execution not in _EXECUTIONS:
            raise ValueError("execution must be 'auto', 'full', 'cutoff', or 'cover'")
        if radius is not None:
            if not isinstance(radius, numbers.Real) or isinstance(radius, (bool, np.bool)):
                raise TypeError("radius must be a real number or None")
            radius = float(radius)

        normalized_cutoff = (
            None
            if cutoff_threshold is None
            else _normalized_nonnegative_integer(cutoff_threshold, "cutoff_threshold")
        )
        normalized_cover = (
            None
            if cover_threshold is None
            else _normalized_nonnegative_integer(cover_threshold, "cover_threshold")
        )
        normalized_threads = _normalized_nonnegative_integer(threads, "threads")
        normalized_fixed_ions: tuple[str, ...]
        if fixed_ions is None:
            normalized_fixed_ions = ()
        else:
            if isinstance(fixed_ions, (str, bytes)):
                raise TypeError("fixed_ions must be a sequence of strings")
            try:
                normalized_fixed_ions = tuple(fixed_ions)
            except TypeError as error:
                raise TypeError("fixed_ions must be a sequence of strings") from error
            if any(not isinstance(ion, str) for ion in normalized_fixed_ions):
                raise TypeError("fixed_ions must contain only strings")

        object.__setattr__(self, "method", method_id)
        object.__setattr__(self, "parameter_set", parameter_set_id)
        object.__setattr__(
            self,
            "options_by_method",
            _frozen_options_by_method(normalized_options),
        )
        object.__setattr__(self, "parameter_matching", parameter_matching)
        object.__setattr__(self, "execution", execution)
        object.__setattr__(self, "radius", radius)
        object.__setattr__(self, "cutoff_threshold", normalized_cutoff)
        object.__setattr__(self, "cover_threshold", normalized_cover)
        object.__setattr__(self, "threads", normalized_threads)
        object.__setattr__(self, "fixed_ions", normalized_fixed_ions)

    @property
    def _permissive_types(self) -> bool:
        return self.parameter_matching == "permissive"
