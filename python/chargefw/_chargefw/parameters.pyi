from collections.abc import Sequence

from .._payloads import ParameterSetDescriptorPayload
from .._types import AtomClassifier, BondClassifier
from .core import _NativeMolecule

class _NativeParameterCatalog:
    def _descriptors(self) -> list[ParameterSetDescriptorPayload]: ...

def _load_bundled_parameter_catalog() -> _NativeParameterCatalog: ...
def _load_parameter_set_catalog(path: str) -> _NativeParameterCatalog: ...
def _make_parameter_template_catalog(
    molecules: Sequence[_NativeMolecule],
    method_id: str,
    atom_types: AtomClassifier,
    bond_types: BondClassifier,
    id: str,
) -> _NativeParameterCatalog: ...
