"""Following a running fit, keeping its best values, and stopping it."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from numpy.typing import NDArray

from .._parameters import ParameterSet
from ._template import _with_values


@dataclass(frozen=True, slots=True, eq=False)
class FitProgress:
    """A snapshot after one trial point of a running fit.

    ``evaluation`` counts the trial points so far, failed ones included; ``failed_evaluations``
    counts those whose calculation failed. ``objective`` is the weighted sum of squared residuals the
    fit minimizes, and ``rmsd`` the RMSD without weights of the charges that have a reference, both
    at the trial point; ``nan`` when its calculation failed. ``best_values`` holds every slot in
    layout order at the lowest objective so far, and ``best_objective`` and ``best_rmsd`` the
    objective and RMSD there; ``best_rmsd`` need not be the lowest RMSD seen. Until a trial point
    succeeds, ``best_values`` holds the start and ``best_objective`` is infinite. ``elapsed_seconds``
    counts from the call of ``fit()``. ``parameter_set()`` turns the best values into a parameter
    set, to save as a checkpoint or to start another fit from.
    """

    evaluation: int
    failed_evaluations: int
    objective: float
    rmsd: float
    best_objective: float
    best_rmsd: float
    best_values: NDArray[np.float64]
    elapsed_seconds: float
    _parameter_set: ParameterSet = field(repr=False, compare=False)

    def parameter_set(
        self, *, id: str, name: str = "", publication: str = "", notes: str = ""
    ) -> ParameterSet:
        """The best values so far as a parameter set for the evaluator's method."""

        return _with_values(self._parameter_set, self.best_values, id, name, publication, notes)


class FitObserver:
    """Receives a ``FitProgress`` after every trial point of a fit and may stop it.

    ``fit()`` calls ``on_evaluation()`` once per trial point, then ``cancelled()``. When that returns
    true, the fit stops; ``fit()`` describes what it returns then.
    """

    def on_evaluation(self, progress: FitProgress) -> None:
        """Handle the snapshot of one evaluation."""

    def cancelled(self) -> bool:
        """Return whether the fit should stop."""

        return False
