"""Default fitting settings for each method: bounds, starting values, and fixed parameters.

``fit()`` uses the entry in ``PROFILES`` for the method being fitted, for whatever the call does not
give. Change a setting with ``dataclasses.replace()`` on a profile, and pass the result as
``fit(..., profile=...)``.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass, field

__all__ = ["PROFILES", "Profile", "sqe_bounds"]


@dataclass(frozen=True, slots=True)
class Profile:
    """What a fit uses for one method unless told otherwise.

    ``bounds`` maps a parameter name to its lower and upper value, for every slot of that name.
    ``start`` maps a parameter name to a starting value, used for a parameter whose values are all
    zero in the plan's parameter set, as in a template. ``fixed`` names parameters whose first slot
    is held at its start value: one per invariance of the method, a change of values that changes no
    charge.
    """

    bounds: Mapping[str, tuple[float, float]] = field(default_factory=dict)
    start: Mapping[str, float] = field(default_factory=dict)
    fixed: tuple[str, ...] = ()


def sqe_bounds(min_width: float = 0.5, min_kappa: float = 0.01) -> dict[str, tuple[float, float]]:
    """Bounds that keep the SQE system positive definite for every molecule and geometry.

    A width of at least ``min_width``, a hardness of at least ``1 / (sqrt(pi) * min_width)``, and a
    bond ``kappa`` of at least ``min_kappa`` keep the smallest eigenvalue of the SQE system at
    ``min_kappa`` or above. A larger ``min_kappa`` keeps the system further from singular but holds
    back charge transfer along bonds, so it is a choice for each fit.
    """

    if not min_width > 0.0 or not min_kappa > 0.0:
        raise ValueError("min_width and min_kappa must be positive")
    return {
        "width": (min_width, math.inf),
        "hardness": (1.0 / (math.sqrt(math.pi) * min_width), math.inf),
        "kappa": (min_kappa, math.inf),
    }


_SQE_START = {"hardness": 2.0, "width": 1.0, "kappa": 1.0}

# The profiles fit() uses by default, by method ID. The fixed slots remove the invariances: in EEM,
# adding a constant to every A and scaling A, B, and kappa together; in the SQE family, adding a
# constant to every electronegativity, and in SQE+qp also to every q0.
PROFILES: Mapping[str, Profile] = {
    # Published EEM sets have a median B / kappa of about 1.5 to 2.2.
    "eem": Profile(
        bounds={"B": (1e-3, math.inf)},
        start={"kappa": 1.0, "B": 2.0},
        fixed=("kappa", "A"),
    ),
    "sqe": Profile(bounds=sqe_bounds(), start=_SQE_START, fixed=("electronegativity",)),
    "sqeq0": Profile(bounds=sqe_bounds(), start=_SQE_START, fixed=("electronegativity",)),
    "sqeqp": Profile(
        bounds=sqe_bounds(),
        start={**_SQE_START, "q0": 0.0},
        fixed=("electronegativity", "q0"),
    ),
}
