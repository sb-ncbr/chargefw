"""The fit's least-squares model from its normal equations (chargefw.parameterisation._compression)."""

from collections.abc import Callable

import numpy as np
import pytest
from chargefw.parameterisation._compression import CompressionError, compress
from numpy.typing import NDArray
from scipy.optimize import least_squares  # type: ignore[import-untyped, unused-ignore]

Vector = NDArray[np.float64]

SETTINGS = {"method": "trf", "x_scale": "jac", "loss": "linear", "tr_solver": "exact"}


def steps(
    residuals: Callable[[Vector], Vector],
    jacobian: Callable[[Vector], Vector],
    start: Vector,
    compressed: bool,
) -> tuple[list[Vector], Vector]:
    """The points SciPy's TRF evaluates, and where it stops, on the residuals themselves or on their
    compressed form."""

    points: list[Vector] = []

    def represented(x: Vector) -> tuple[Vector, Vector]:
        r, J = residuals(x), jacobian(x)
        return compress(float(r @ r), J.T @ J, J.T @ r)

    def fun(x: Vector) -> Vector:
        points.append(x.copy())
        return represented(x)[0] if compressed else residuals(x)

    def jac(x: Vector) -> Vector:
        return represented(x)[1] if compressed else jacobian(x)

    result = least_squares(fun, start, jac=jac, **SETTINGS)
    return points, result.x


@pytest.mark.parametrize(
    ("residuals", "jacobian", "start"),
    [
        # A slot with a tiny derivative that still matters: the minimum is at x2 = -1e7.
        (
            lambda x: np.array([x[0], 1.0 + 1e-7 * x[1]]),
            lambda x: np.array([[1.0, 0.0], [0.0, 1e-7]]),
            np.zeros(2),
        ),
        # Rank 1 at the start, full rank afterwards, and a residual that reaches zero.
        (
            lambda x: np.array([x[0] * x[1] - 1.0, x[1] - 0.5]),
            lambda x: np.array([[x[1], x[0]], [0.0, 1.0]]),
            np.array([1.0, 0.0]),
        ),
        # Slots in very different units.
        (
            lambda x: np.array([10.0 * (x[1] / 1e3 - (x[0] / 1e-3) ** 2), 1.0 - x[0] / 1e-3]),
            lambda x: np.array([[-20.0 * x[0] / 1e-6, 10.0 / 1e3], [-1.0 / 1e-3, 0.0]]),
            np.array([-1.2e-3, 1e3]),
        ),
    ],
)
def test_compressed_model_takes_the_steps_of_the_residuals(
    residuals: Callable[[Vector], Vector], jacobian: Callable[[Vector], Vector], start: Vector
) -> None:
    dense_points, dense_x = steps(residuals, jacobian, start, compressed=False)
    compressed_points, compressed_x = steps(residuals, jacobian, start, compressed=True)

    assert len(compressed_points) == len(dense_points)
    for dense, compressed in zip(dense_points, compressed_points, strict=True):
        np.testing.assert_allclose(compressed, dense, rtol=1e-9, atol=1e-12)
    np.testing.assert_allclose(compressed_x, dense_x, rtol=1e-9)


def test_compressed_model_reproduces_the_normal_equations() -> None:
    rng = np.random.default_rng(1)
    J = rng.normal(size=(40, 4)) * np.array([1.0, 1e-5, 1e3, 1.0])
    J[:, 3] = J[:, 0] + J[:, 1]  # a null direction
    r = rng.normal(size=40)
    r -= J @ np.linalg.lstsq(J, r, rcond=None)[0] * 0.5
    F, H, z = float(r @ r), J.T @ J, J.T @ r

    f, matrix = compress(F, H, z)

    assert f.shape == (5,) and matrix.shape == (5, 4)
    np.testing.assert_allclose(matrix.T @ matrix, H, rtol=1e-10, atol=1e-10 * np.abs(H).max())
    np.testing.assert_allclose(matrix.T @ f, z, rtol=1e-8, atol=1e-8 * np.abs(z).max())
    assert f @ f == pytest.approx(F, rel=1e-12)


@pytest.mark.parametrize(
    ("total", "gram", "gradient"),
    [
        (1.0, np.array([[1.0, 2.0], [2.0, 1.0]]), np.array([0.1, 0.1])),  # indefinite
        (0.4, np.eye(2), np.array([1.0, 0.0])),  # the gradient implies |f|^2 = 1, not 0.4
        # Nearly dependent columns: the weakest direction is kept, but its eigenvector is too
        # inaccurate for its part of the gradient, which would overstate F = 1.
        (
            1.0,
            np.array([[1.0, 1.0], [1.0, 1.0 + 3.8062979230107953e-7**2]]),
            np.array([0.0, 3.8062979230107953e-7]),
        ),
        (1.0, np.diag([1.0, 0.0]), np.array([0.1, 0.2])),  # a gradient without derivatives
        (np.nan, np.eye(2), np.zeros(2)),
        (1.0, np.diag([1.0, -1e-3]), np.zeros(2)),
    ],
)
def test_compression_rejects_inconsistent_normal_equations(
    total: float, gram: Vector, gradient: Vector
) -> None:
    with pytest.raises(CompressionError):
        compress(total, gram, gradient)
