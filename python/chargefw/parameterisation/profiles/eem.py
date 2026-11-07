"""The fitting profile of EEM."""

from __future__ import annotations

import math

from ._profile import Profile

__all__ = ["PROFILES"]

# Published EEM sets have a median B / kappa of about 1.5 to 2.2. The fixed slots remove the
# invariances: a constant added to every A, and scaling A, B, and kappa together.
PROFILES = {
    "eem": Profile(
        bounds={"B": (1e-3, math.inf)},
        start={"kappa": 1.0, "B": 2.0},
        fixed=("kappa", "A"),
    ),
}
