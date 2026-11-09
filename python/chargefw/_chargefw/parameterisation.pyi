from collections.abc import Sequence

import numpy as np
from numpy.typing import NDArray

from .._types import AtomClassifier, BondClassifier, ParameterScope
from ..parameterisation._payloads import CandidateChargesPayload, NormalEquationsPayload
from .calculation import _NativePlan
from .core import _NativeMolecule
from .parameters import _NativeParameterCatalog

class _NativeEvaluator:
    def __init__(self, plan: _NativePlan) -> None: ...
    def _layout(self) -> list[tuple[ParameterScope, int, str]]: ...
    def _values(self) -> NDArray[np.float64]: ...
    def _calculate_values(
        self, values: NDArray[np.float64], max_threads: int | None = None
    ) -> CandidateChargesPayload: ...
    def _entry_counts(self) -> tuple[list[int], list[int]]: ...
    def _single_bond_entries(self) -> list[tuple[bool, bool, list[int]]]: ...
    def _normal_equations(
        self,
        values: NDArray[np.float64],
        references: NDArray[np.float64],
        weights: NDArray[np.float64],
        max_threads: int | None = None,
    ) -> NormalEquationsPayload: ...
    def _exact_derivatives(self) -> bool: ...

def _make_parameter_template_catalog(
    molecules: Sequence[_NativeMolecule],
    method_id: str,
    atom_types: AtomClassifier,
    bond_types: BondClassifier,
    id: str,
) -> _NativeParameterCatalog: ...
def _with_values(
    catalog: _NativeParameterCatalog,
    id: str,
    values: NDArray[np.float64],
    new_id: str,
    name: str,
    publication: str,
    notes: str,
) -> _NativeParameterCatalog: ...
