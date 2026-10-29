from typing import Literal, TypedDict

from .calculation import _NativeExecutionResult

class _NativeInputMetadata: ...

class HierarchyLabelsPayload(TypedDict):
    atom: str | None
    residue: str | None
    chain: str | None
    sequence: str | None

class StructuralLabelsPayload(TypedDict):
    author: HierarchyLabelsPayload
    label: HierarchyLabelsPayload
    entity: str | None
    insertion_code: str | None
    alternate_location: str | None
    segment: str | None

type AtomReferencePayload = tuple[int, str | None, StructuralLabelsPayload | None]
type ConformerReferencePayload = tuple[int, str | None, list[AtomReferencePayload]]
type ComponentInstancePayload = tuple[str, list[int]]

class ImportMetadataPayload(TypedDict):
    format: Literal["mol", "sdf", "mol2", "molecule-json", "pdb", "mmcif"]
    atoms: list[AtomReferencePayload]
    conformers: list[ConformerReferencePayload]
    record_selection: str | None
    alternate_location_selection: str | None
    conformer_selection: str | None
    bond_strategy: str | None
    source_connectivity: Literal["absent", "explicitly-empty", "present"]
    components: list[ComponentInstancePayload]

class MoleculePayload(TypedDict):
    atomic_numbers: list[int]
    formal_charges: list[int]
    bonds: list[tuple[int, int, int]]
    coordinates: list[list[tuple[float, float, float]]]
    name: str
    atom_names: list[str]
    conformer_names: list[str]
    source: str
    record_index: int
    record_id: str | int | None
    import_metadata: ImportMetadataPayload | None
    native_input_metadata: _NativeInputMetadata

def _parse(
    contents: str,
    source: str,
    format: Literal["mol", "sdf", "mol2", "molecule-json", "pdb", "mmcif"],
    selection: Literal["all", "polymers-and-ligands", "polymers"],
    bonds: Literal["none", "explicit", "templates", "hybrid"],
    conformers: Literal["first", "all"],
) -> list[MoleculePayload]: ...
def _dumps(
    results: list[_NativeExecutionResult],
    format: Literal["mol2", "mmcif", "pqr", "result-json"],
) -> str: ...
def _attach_mmcif(contents: str, results: list[_NativeExecutionResult], overwrite: bool) -> str: ...
