"""The fit's least-squares model in the size of the free slots instead of the number of charges.

From the normal equations of the weighted residuals at one point, F = |r|², H = JᵀJ and z = Jᵀr in
n free slots, compress() builds n + 1 residuals f and their (n + 1) × n matrix J̃ with |f|² = F,
and with J̃ᵀJ̃ = H and J̃ᵀf = z except in the dropped directions described below. SciPy's TRF method
with the linear loss uses only these quantities, so it sees the local quadratic model of the full
problem in the directions kept; with directions dropped, its steps can differ from those on the full
problem. J̃ is not the derivative of f; it represents the model at that point only.

H is scaled to unit diagonal first, so that a slot with small derivatives keeps its direction. Its
directions with an eigenvalue at or below undetermined() of the largest are dropped: combinations of
slots that the data determines less than about 1e-6 as strongly as the best determined one, in the
coordinates the optimizer uses. The fit's diagnostics count such combinations by the same threshold,
but in the parameter values themselves, so the two usually, not always, agree. The optimizer takes no
step to reduce the
residual along them, which differs from plain least squares: that would follow such a combination,
usually to extreme parameter values. Its steps can still move along them within bounds.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

__all__ = ["CompressionError", "compress", "determined_step", "undetermined"]

_EPSILON = float(np.finfo(np.float64).eps)
# How much more than F the residuals may explain, relative to F, at least: rounding, not a different
# model.
_SUM_OF_SQUARES_TOLERANCE = 1e-8


class CompressionError(ArithmeticError):
    """The normal equations cannot represent the least-squares model in floating point."""


def undetermined(size: int) -> float:
    """The relative eigenvalue of a unit-diagonal Gram matrix in `size` slots at or below which a
    direction counts as not determined: singular values below 1e-6 of the largest, or below what
    floating point resolves when that is more, as with thousands of slots."""

    return max(1e-12, 64.0 * size * _EPSILON)


def _decompose(
    gram: NDArray[np.float64],
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], NDArray[np.bool_], float]:
    """H scaled to unit diagonal, S H S, decomposed: the scale S per slot, the eigenvalues and
    eigenvectors, the directions undetermined(), and the threshold for them."""

    diagonal = np.diag(gram)
    if np.any(diagonal < 0.0):
        raise CompressionError("the normal equations have a negative diagonal")
    norms = np.sqrt(diagonal)
    scale = np.ones(diagonal.size)
    np.divide(1.0, norms, out=scale, where=norms != 0.0)
    scaled = scale[:, np.newaxis] * gram * scale[np.newaxis, :]
    eigenvalues, vectors = np.linalg.eigh(scaled)
    threshold = undetermined(diagonal.size) * max(float(eigenvalues[-1]), 0.0)
    if eigenvalues[0] < -threshold:
        raise CompressionError(
            f"the normal equations are not positive semidefinite (eigenvalue {eigenvalues[0]:.3e})"
        )
    return scale, eigenvalues, vectors, eigenvalues <= threshold, threshold


def determined_step(
    gram: NDArray[np.float64], gradient: NDArray[np.float64]
) -> NDArray[np.float64]:
    """The Gauss-Newton step in the directions compress() keeps, and no step in the others: the
    shortest step to the least-squares minimum of the linear model in the coordinates scaled to unit
    diagonal, which is in general not the shortest one in the slots themselves. For residuals that
    are affine in the slots, it reaches their least-squares minimum along the directions the data
    determines; components in the other directions keep their start."""

    if not (np.all(np.isfinite(gram)) and np.all(np.isfinite(gradient))):
        raise CompressionError("the normal equations are not finite")
    scale, eigenvalues, vectors, dropped, _ = _decompose(gram)
    kept = vectors[:, ~dropped]
    # In scaled coordinates y = S⁻¹x the Gram matrix is V Λ Vᵀ; the step there is -V Λ⁻¹ Vᵀ S z on
    # the kept directions, and x = S y.
    step = kept @ ((kept.T @ (scale * gradient)) / eigenvalues[~dropped])
    return -scale * step


def compress(
    sum_of_squares: float, gram: NDArray[np.float64], gradient: NDArray[np.float64]
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    n = gradient.size
    if not (
        np.isfinite(sum_of_squares) and np.all(np.isfinite(gram)) and np.all(np.isfinite(gradient))
    ):
        raise CompressionError("the normal equations are not finite")
    scale, eigenvalues, vectors, dropped, threshold = _decompose(gram)
    if np.any(gradient[np.diag(gram) == 0.0] != 0.0):
        raise CompressionError("the gradient is not zero for a slot without derivatives")
    scaled_gradient = scale * gradient

    # A direction v with singular value σ carries at most σ·√F of the gradient: |vᵀz| = |(Jv)ᵀr|
    # ≤ |Jv|·|r| (Cauchy-Schwarz). Dropped directions have σ² ≤ threshold, so more than that, with
    # room for rounding, means normal equations that do not belong together.
    if np.any(dropped):
        projected = float(np.linalg.norm(vectors[:, dropped].T @ scaled_gradient))
        bound = np.sqrt(float(np.sum(dropped)) * threshold * sum_of_squares)
        if projected > 4.0 * bound + 64.0 * _EPSILON * float(np.linalg.norm(scaled_gradient)):
            raise CompressionError(
                "the gradient in directions the normal equations cannot resolve exceeds what their "
                f"derivatives allow ({projected:.3e} > {bound:.3e})"
            )

    root = np.sqrt(np.where(dropped, 0.0, eigenvalues))
    small = np.where(dropped, 0.0, (vectors.T @ scaled_gradient) / np.where(dropped, 1.0, root))
    # R = Λ^(1/2) Vᵀ S⁻¹, so RᵀR = S⁻¹ V Λ Vᵀ S⁻¹ = H on the directions kept.
    matrix = root[:, np.newaxis] * vectors.T
    np.divide(matrix, scale[np.newaxis, :], out=matrix)
    # In exact arithmetic |small|² is the part of F the model explains, at most F. Rounding in the
    # eigenvectors of the kept directions, about n·ε times the largest over the smallest kept
    # eigenvalue, can push it above F by that much, and the gradient Rᵀ·small is then off by as much.
    # Within that, small is scaled down to |small|² = F, so that the optimizer sees the exact sum of
    # squares and a gradient within its rounding. More means normal equations that do not belong
    # together.
    explained = float(small @ small)
    kept = eigenvalues[~dropped]
    rounding = 64.0 * n * _EPSILON * float(eigenvalues[-1] / kept[0]) if kept.size else 0.0
    if explained > sum_of_squares * (1.0 + max(_SUM_OF_SQUARES_TOLERANCE, rounding)):
        raise CompressionError(
            f"the gradient implies a sum of squares of {explained:.6e}, not {sum_of_squares:.6e}"
        )
    if explained > sum_of_squares:
        small *= np.sqrt(sum_of_squares / explained)
        explained = float(small @ small)
    remainder = sum_of_squares - explained
    residuals = np.append(small, np.sqrt(max(remainder, 0.0)))
    return residuals, np.vstack([matrix, np.zeros(n)])
