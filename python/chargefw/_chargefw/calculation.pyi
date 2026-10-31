from collections.abc import Sequence

from .._payloads import AssessmentReportPayload, ExecutionPlanPayload, ExecutionResultPayload
from .._types import PortableId, RequestedExecution
from .adapters import _NativeInputMetadata
from .core import _NativeMolecule
from .parameters import _NativeParameterCatalog

def _fixed_ion_names(common_only: bool) -> tuple[str, ...]: ...

class _NativeExecutionResult:
    def report(self) -> ExecutionResultPayload: ...

class _NativePlan:
    def report(self) -> ExecutionPlanPayload: ...
    def calculate(
        self, max_threads: int | None = None, observer: object | None = None
    ) -> _NativeExecutionResult: ...
    def estimate(self, max_threads: int | None = None) -> tuple[float, float]: ...

class _NativeAssessment:
    def report(self) -> AssessmentReportPayload: ...
    def plans(self) -> list[_NativePlan]: ...
    def calculate_default(self, observer: object | None = None) -> _NativeExecutionResult: ...

def _make_assessment(
    molecules: Sequence[_NativeMolecule],
    input_metadata: Sequence[_NativeInputMetadata | None],
    identities: Sequence[tuple[str, int, str | int | None]],
    caller_atom_ids: Sequence[Sequence[PortableId] | None],
    molecule_collection_name: str,
    catalog: _NativeParameterCatalog,
    method_id: str | None,
    parameter_set_id: str | None,
    method_options: dict[str, dict[str, bool | int | float | str]],
    permissive_types: bool,
    execution: RequestedExecution,
    radius: float | None,
    cutoff_threshold: int | None,
    cover_threshold: int | None,
    max_threads: int,
    fixed_ions: Sequence[str] = (),
) -> _NativeAssessment: ...
