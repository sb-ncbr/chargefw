"""What a fit's result shows about how far its values can be trusted."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from types import MappingProxyType

import numpy as np
from numpy.typing import NDArray

from .._elements import periodic_table
from ._compression import undetermined
from ._evaluator import EvaluationError


@dataclass(frozen=True, slots=True)
class ElementDiagnostics:
    """How well the charges of one element match the references.

    ``r2`` is ``nan`` when every reference of the element is the same, which leaves it undefined.
    """

    count: int
    rmsd: float
    r2: float


@dataclass(frozen=True, slots=True)
class Diagnostics:
    """What shows whether the fitted values can be trusted.

    ``r2`` is the coefficient of determination of the charges at the fitted values against the
    references, without weights; ``nan`` when every reference is the same. ``at_bounds`` lists free
    slots whose fitted value is at a bound, within ``1e-8 * max(1, |bound|)``; for a slot the profile
    transforms, the bound of the transformed quantity counts, such as the stability limit of the SQE
    hardness. ``failed_messages`` holds the distinct messages of trial calculations that failed, at
    most 20.

    ``singular_values`` come from the weighted Jacobian at the fitted values, for the free slots,
    with each column scaled to unit length, so they compare slots whatever their units; they are in
    decreasing order. ``rank`` counts those above ``1e-6`` times the largest (somewhat more with
    thousands of free slots, where rounding is larger), well above the noise of a finite-difference
    Jacobian. A rank below the number of free slots means the data cannot determine some combination
    of them; a fit with exact derivatives leaves out combinations it finds this weak. ``weakest``
    lists the free slots that dominate the direction of the smallest singular value: those whose
    share of it is at least half the largest share.

    ``matches`` gives, for each free slot in the order of ``FitResult.free``, the number of
    calculated atoms or bonds of its type, or ``None`` for a common parameter. ``by_element`` maps
    each element symbol, in atomic-number order, to the number of charges of that element and their
    RMSD and R² against the references.
    """

    r2: float
    by_element: Mapping[str, ElementDiagnostics]
    matches: tuple[int | None, ...]
    at_bounds: tuple[int, ...]
    failed_messages: tuple[str, ...]
    singular_values: NDArray[np.float64]
    rank: int
    weakest: tuple[int, ...]


def _spectrum_of_jacobian(
    jacobian: NDArray[np.float64],
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """The singular values of the Jacobian with unit columns, in decreasing order, one per column,
    and the right singular vector of the smallest."""

    norms = np.linalg.norm(jacobian, axis=0)
    scaled = jacobian / np.where(norms > 0.0, norms, 1.0)
    # The right singular vectors give the weakest direction. A reduced SVD has all of them when there
    # are at least as many charges as free slots; only an underdetermined Jacobian needs the full one
    # for its null directions, and then the left vectors are small. A full SVD of a tall Jacobian
    # would build a square matrix in the number of charges.
    rows, columns = scaled.shape
    _, singular_values, right = np.linalg.svd(scaled, full_matrices=rows < columns)
    singular_values = np.concatenate((singular_values, np.zeros(columns - singular_values.size)))
    return singular_values, right[-1]


def _spectrum_of_gram(
    gram: NDArray[np.float64],
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """As _spectrum_of_jacobian(), from JᵀJ: its eigenvalues, scaled to unit diagonal, are the
    squared singular values, and its eigenvectors the right singular vectors."""

    norms = np.sqrt(np.clip(np.diag(gram), 0.0, None))
    scale = np.ones(norms.size)
    np.divide(1.0, norms, out=scale, where=norms > 0.0)
    eigenvalues, vectors = np.linalg.eigh(scale[:, np.newaxis] * gram * scale[np.newaxis, :])
    if eigenvalues[0] < -undetermined(norms.size) * max(float(eigenvalues[-1]), 0.0):
        raise EvaluationError(
            "the normal equations at the fitted values are not positive semidefinite "
            f"(eigenvalue {eigenvalues[0]:.3e})"
        )
    return np.sqrt(np.clip(eigenvalues, 0.0, None))[::-1].copy(), vectors[:, 0]


def _diagnostics(
    charges: NDArray[np.float64],
    reference: NDArray[np.float64],
    spectrum: tuple[NDArray[np.float64], NDArray[np.float64]],
    values: NDArray[np.float64],
    free: Sequence[int],
    lower: NDArray[np.float64],
    upper: NDArray[np.float64],
    failed_messages: Sequence[str],
    elements: NDArray[np.int64],
    matches: tuple[int | None, ...],
) -> Diagnostics:
    by_element = {
        periodic_table[int(atomic_number)].symbol: ElementDiagnostics(
            count=int(np.count_nonzero(elements == atomic_number)),
            rmsd=_rmsd(charges[elements == atomic_number], reference[elements == atomic_number]),
            r2=_r2(charges[elements == atomic_number], reference[elements == atomic_number]),
        )
        for atomic_number in np.unique(elements)
    }
    r2 = _r2(charges, reference)

    def at_bound(value: float, bound: float) -> bool:
        return math.isfinite(bound) and abs(value - bound) <= 1e-8 * max(1.0, abs(bound))

    at_bounds = tuple(
        slot
        for slot in free
        if at_bound(values[slot], lower[slot]) or at_bound(values[slot], upper[slot])
    )

    singular_values, smallest_direction = spectrum
    # The threshold by which the compressed fit leaves combinations out, here applied to the
    # parameter values rather than the optimizer's coordinates.
    threshold = np.sqrt(undetermined(singular_values.size)) * singular_values[0]
    rank = int(np.sum(singular_values > threshold)) if singular_values[0] else 0
    weakest_direction = np.abs(smallest_direction)
    weakest = tuple(
        free[index]
        for index in np.argsort(-weakest_direction, kind="stable")
        if weakest_direction[index] >= 0.5 * weakest_direction.max()
    )
    singular_values.setflags(write=False)
    return Diagnostics(
        r2=r2,
        by_element=MappingProxyType(by_element),
        matches=matches,
        at_bounds=at_bounds,
        failed_messages=tuple(failed_messages),
        singular_values=singular_values,
        rank=rank,
        weakest=weakest,
    )


def _rmsd(charges: NDArray[np.float64], reference: NDArray[np.float64]) -> float:
    return float(np.sqrt(np.mean((charges - reference) ** 2)))


def _r2(charges: NDArray[np.float64], reference: NDArray[np.float64]) -> float:
    total = float(np.sum((reference - reference.mean()) ** 2))
    residual = float(np.sum((charges - reference) ** 2))
    return 1.0 - residual / total if total > 0.0 else math.nan
