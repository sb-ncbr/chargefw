"""Private normalization for Python-facing calculation arguments."""

from __future__ import annotations

import numbers
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from types import MappingProxyType
from typing import SupportsIndex

import numpy as np

from ._arrays import as_integer
from ._chargefw import calculation as _native_calculation
from ._methods import Method
from ._parameters import ParameterSet
from ._types import MethodOptionValue, MethodScope, ParameterMatching, RequestedExecution

_EXECUTIONS = frozenset(("auto", "full", "cutoff", "cover"))


def _normalized_nonnegative_integer(value: SupportsIndex, name: str) -> int:
    result = as_integer(value, name)
    if result < 0:
        raise ValueError(f"{name} must be non-negative")
    if result > np.iinfo(np.uintp).max:
        raise ValueError(f"{name} is outside the native size range")
    return result


@dataclass(frozen=True, slots=True, init=False)
class FixedAtomCharge:
    """A fixed charge on one atom; component_id is the residue name it was selected by."""

    molecule_index: int
    atom_index: int
    charge: float
    component_id: str | None

    def __init__(
        self,
        molecule_index: SupportsIndex,
        atom_index: SupportsIndex,
        charge: float,
        component_id: str | None = None,
    ) -> None:
        if isinstance(charge, (bool, np.bool)) or not isinstance(charge, numbers.Real):
            raise TypeError("charge must be a real number")
        if component_id is not None and not isinstance(component_id, str):
            raise TypeError("component_id must be a string or None")
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
        object.__setattr__(self, "component_id", component_id)


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


@dataclass(frozen=True, slots=True, init=False)
class RequestedCalculation:
    """Normalized policy requested for an assessment or calculation."""

    method: Method | None
    parameter_set: ParameterSet | None
    method_scope: MethodScope
    method_options: Mapping[str, MethodOptionValue]
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
        method: Method | None = None,
        parameter_set: ParameterSet | None = None,
        method_scope: MethodScope = "recommended",
        method_options: Mapping[str, MethodOptionValue] | None = None,
        parameter_matching: ParameterMatching = "strict",
        execution: RequestedExecution = "auto",
        radius: float | None = None,
        cutoff_threshold: int | None = _native_calculation._default_cutoff_atom_threshold,
        cover_threshold: int | None = _native_calculation._default_cover_atom_threshold,
        threads: int = 0,
        fixed_ions: Sequence[str] | None = None,
    ) -> None:
        if method is not None and not isinstance(method, Method):
            raise TypeError("method must be a Method or None")
        if parameter_set is not None and not isinstance(parameter_set, ParameterSet):
            raise TypeError("parameter_set must be a ParameterSet or None")
        method_id = None if method is None else method.id
        if not isinstance(method_scope, str):
            raise TypeError("method_scope must be a string")
        if method_scope not in ("recommended", "all"):
            raise ValueError("method_scope must be 'recommended' or 'all'")
        if method_scope == "all" and method_id is not None:
            raise ValueError("method_scope applies only without an explicit method")
        # Options belong to the named method: with automatic selection, an option such as
        # iters would apply to whichever method happened to be chosen.
        normalized_options: dict[str, MethodOptionValue] = {}
        if method_options is not None:
            if not isinstance(method_options, Mapping):
                raise TypeError("method_options must be a mapping")
            if method_options:
                if method_id is None:
                    raise ValueError("method_options requires an explicit method")
                normalized_options = _normalized_option_values(method_options, method_id)

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

        object.__setattr__(self, "method", method)
        object.__setattr__(self, "parameter_set", parameter_set)
        object.__setattr__(self, "method_scope", method_scope)
        object.__setattr__(self, "method_options", MappingProxyType(normalized_options))
        object.__setattr__(self, "parameter_matching", parameter_matching)
        object.__setattr__(self, "execution", execution)
        object.__setattr__(self, "radius", radius)
        object.__setattr__(self, "cutoff_threshold", normalized_cutoff)
        object.__setattr__(self, "cover_threshold", normalized_cover)
        object.__setattr__(self, "threads", normalized_threads)
        object.__setattr__(self, "fixed_ions", normalized_fixed_ions)
