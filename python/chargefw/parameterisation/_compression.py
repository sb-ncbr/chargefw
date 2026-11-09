"""The fit's least-squares model in the size of the free slots instead of the number of charges.

From the normal equations of the weighted residuals at one point, F = |r|², H = JᵀJ and z = Jᵀr in
n free slots, compress() builds n + 1 residuals f and their (n + 1) × n matrix J̃ with J̃ᵀJ̃ = H,
J̃ᵀf = z and |f|² = F. SciPy's TRF method with the linear loss uses only these quantities, so it
takes the steps of the full problem. J̃ is not the derivative of f; it represents the model at that
point only.

H is scaled to unit diagonal first, so that a slot with small derivatives keeps its direction. Only
directions that the scaled H cannot resolve in floating point are dropped, and only when the
gradient has no material part in them.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray

__all__ = ["CompressionError", "compress"]

_EPSILON = float(np.finfo(np.float64).eps)
# The largest share of the scaled gradient allowed in directions that are dropped.
_GRADIENT_TOLERANCE = 1e-6
# How much more than F the residuals may explain, relative to F: rounding, not a different model.
_SUM_OF_SQUARES_TOLERANCE = 1e-8


class CompressionError(ArithmeticError):
    """The normal equations cannot represent the least-squares model in floating point."""


def resolution(size: int) -> float:
    """Eigenvalues of a unit-diagonal Gram matrix below this times the largest are rounding."""

    return 64.0 * size * _EPSILON


def compress(
    sum_of_squares: float, gram: NDArray[np.float64], gradient: NDArray[np.float64]
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    n = gradient.size
    if not (
        np.isfinite(sum_of_squares) and np.all(np.isfinite(gram)) and np.all(np.isfinite(gradient))
    ):
        raise CompressionError("the normal equations are not finite")
    diagonal = np.diag(gram)
    if np.any(diagonal < 0.0):
        raise CompressionError("the normal equations have a negative diagonal")
    norms = np.sqrt(diagonal)
    without = norms == 0.0
    if np.any(gradient[without] != 0.0):
        raise CompressionError("the gradient is not zero for a slot without derivatives")
    scale = np.ones(n)
    np.divide(1.0, norms, out=scale, where=~without)
    scaled = scale[:, np.newaxis] * gram * scale[np.newaxis, :]
    scaled_gradient = scale * gradient

    eigenvalues, vectors = np.linalg.eigh(scaled)
    threshold = resolution(n) * max(float(eigenvalues[-1]), 0.0)
    if eigenvalues[0] < -threshold:
        raise CompressionError(
            f"the normal equations are not positive semidefinite (eigenvalue {eigenvalues[0]:.3e})"
        )
    dropped = eigenvalues <= threshold
    gradient_norm = float(np.linalg.norm(scaled_gradient))
    if np.any(dropped) and gradient_norm > 0.0:
        share = float(np.linalg.norm(vectors[:, dropped].T @ scaled_gradient)) / gradient_norm
        if share > _GRADIENT_TOLERANCE:
            raise CompressionError(
                f"{share:.1e} of the gradient lies in directions the normal equations cannot resolve"
            )

    root = np.sqrt(np.where(dropped, 0.0, eigenvalues))
    small = np.where(dropped, 0.0, (vectors.T @ scaled_gradient) / np.where(dropped, 1.0, root))
    # R = Λ^(1/2) Vᵀ S⁻¹, so RᵀR = S⁻¹ V Λ Vᵀ S⁻¹ = H on the directions kept.
    matrix = root[:, np.newaxis] * vectors.T
    np.divide(matrix, scale[np.newaxis, :], out=matrix)
    # In exact arithmetic |small|² is the part of F the model explains, at most F. More than that
    # means a direction just above the resolution whose eigenvector is too inaccurate for its part of
    # the gradient; neither the remainder nor f can then be adjusted without changing the sum of
    # squares or the gradient the optimizer sees.
    explained = float(small @ small)
    if explained > sum_of_squares * (1.0 + _SUM_OF_SQUARES_TOLERANCE):
        raise CompressionError(
            f"the gradient implies a sum of squares of {explained:.6e}, not {sum_of_squares:.6e}"
        )
    remainder = sum_of_squares - explained
    residuals = np.append(small, np.sqrt(max(remainder, 0.0)))
    return residuals, np.vstack([matrix, np.zeros(n)])
