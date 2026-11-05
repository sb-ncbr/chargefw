from collections.abc import Sequence

import numpy as np
from numpy.typing import NDArray

from .._payloads import ParameterSetDescriptorPayload
from .._types import AtomClassifier, BondClassifier
from .core import _NativeMolecule

class _NativeParameterCatalog:
    def _descriptors(self) -> list[ParameterSetDescriptorPayload]: ...
    def _json(self, id: str) -> str: ...
    def _with_values(
        self,
        id: str,
        values: NDArray[np.float64],
        new_id: str,
        name: str,
        publication: str,
        notes: str,
    ) -> _NativeParameterCatalog: ...

def _load_bundled_parameter_catalog() -> _NativeParameterCatalog: ...
def _load_parameter_set_catalog(path: str) -> _NativeParameterCatalog: ...
def _make_parameter_template_catalog(
    molecules: Sequence[_NativeMolecule],
    method_id: str,
    atom_types: AtomClassifier,
    bond_types: BondClassifier,
    id: str,
) -> _NativeParameterCatalog: ...
