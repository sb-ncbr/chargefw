"""Results of the private native fitting functions (_chargefw.parameterization)."""

from __future__ import annotations

from typing import TypedDict

import numpy as np
from numpy.typing import NDArray

from .._types import ExecutionStatus


class CandidateChargesPayload(TypedDict):
    status: ExecutionStatus
    charges: NDArray[np.float64] | None
    assignment_sizes: list[int] | None
    assignment_molecules: list[int] | None
    failure_message: str | None


class NormalEquationsPayload(TypedDict):
    status: ExecutionStatus
    failure_message: str | None
    charges: NDArray[np.float64] | None
    sum_of_squares: float | None
    # Slots x slots, flattened row by row.
    gauss_newton: NDArray[np.float64] | None
    gradient: NDArray[np.float64] | None
