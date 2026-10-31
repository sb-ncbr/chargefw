"""Method catalog and calculation issue values."""

from __future__ import annotations

from dataclasses import dataclass

from ._catalog import _Catalog
from ._chargefw import methods as _native_methods
from ._parameters import ParameterSetCatalog
from ._payloads import MethodDescriptorPayload
from ._types import (
    ExecutionIssueKind,
    MethodOptionType,
    PrerequisiteIssueKind,
)


@dataclass(frozen=True, slots=True)
class PrerequisiteIssue:
    kind: PrerequisiteIssueKind
    message: str
    molecule_index: int | None = None
    atom_index: int | None = None
    bond_index: int | None = None
    conformer_index: int | None = None


@dataclass(frozen=True, slots=True)
class ExecutionIssue:
    kind: ExecutionIssueKind
    message: str
    molecule_index: int | None = None


@dataclass(frozen=True, slots=True, repr=False)
class MethodOption:
    id: str
    description: str
    type: MethodOptionType
    default: bool | int | float | str
    choices: tuple[bool | int | float | str, ...]
    minimum: bool | int | float | str | None
    minimum_inclusive: bool
    maximum: bool | int | float | str | None
    maximum_inclusive: bool

    def __repr__(self) -> str:
        return f"{type(self).__name__}(id={self.id!r}, default={self.default!r})"


class MethodOptionCatalog(_Catalog[MethodOption]):
    """Immutable ordered method options with lookup by ID."""

    __slots__ = ()

    def __init__(self, values: tuple[MethodOption, ...]) -> None:
        super().__init__(values, "method option")


@dataclass(frozen=True, slots=True, repr=False)
class Method:
    id: str
    name: str
    full_name: str
    publication: str | None
    notes: str
    recommended: bool
    requires_coordinates: bool
    requires_bonds: bool
    time_complexity: str
    memory_complexity: str
    supports_cutoff: bool
    supports_cover: bool
    supports_fixed_ions: bool
    options: MethodOptionCatalog
    parameter_sets: ParameterSetCatalog

    def __repr__(self) -> str:
        return f"{type(self).__name__}(id={self.id!r}, name={self.name!r})"


class MethodCatalog(_Catalog[Method]):
    """Immutable ordered built-in methods with lookup by ID."""

    __slots__ = ()

    def __init__(self, values: tuple[Method, ...]) -> None:
        super().__init__(values, "method")


def _method_descriptor(
    value: MethodDescriptorPayload, parameter_sets: ParameterSetCatalog
) -> Method:
    return Method(
        id=value["id"],
        name=value["name"],
        full_name=value["full_name"],
        publication=value["publication"],
        notes=value["notes"],
        recommended=value["recommended"],
        requires_coordinates=value["requires_coordinates"],
        requires_bonds=value["requires_bonds"],
        time_complexity=value["time_complexity"],
        memory_complexity=value["memory_complexity"],
        supports_cutoff=value["supports_cutoff"],
        supports_cover=value["supports_cover"],
        supports_fixed_ions=value["supports_fixed_ions"],
        options=MethodOptionCatalog(
            tuple(
                MethodOption(
                    id=option["id"],
                    description=option["description"],
                    type=option["type"],
                    default=option["default"],
                    choices=tuple(option["choices"]),
                    minimum=option["minimum"],
                    minimum_inclusive=option["minimum_inclusive"],
                    maximum=option["maximum"],
                    maximum_inclusive=option["maximum_inclusive"],
                )
                for option in value["options"]
            )
        ),
        parameter_sets=ParameterSetCatalog(
            tuple(
                parameter_set
                for parameter_set in parameter_sets.values()
                if parameter_set.method == value["id"]
            )
        ),
    )


def _method_catalog(parameter_sets: ParameterSetCatalog) -> MethodCatalog:
    """Build method descriptors linked to the bundled parameter catalog."""

    return MethodCatalog(
        tuple(
            _method_descriptor(value, parameter_sets)
            for value in _native_methods._method_descriptors()
        )
    )


for _value_type in (
    PrerequisiteIssue,
    ExecutionIssue,
    MethodOption,
    Method,
):
    _value_type.__module__ = "chargefw"
