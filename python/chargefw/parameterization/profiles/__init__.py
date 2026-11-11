"""Default fitting settings for each method: bounds, starting values, fixed parameters, and changes of
variables.

``fit()`` uses the entry in ``PROFILES`` for the method being fitted, for whatever the call does not
give. Change a setting with ``dataclasses.replace()`` on a profile, and pass the result as
``fit(..., profile=...)``. The settings of each method live in its own module (``eem``, ``sqe``).
"""

from __future__ import annotations

from collections.abc import Mapping
from types import MappingProxyType

from . import eem, sqe
from ._profile import Profile, Transform
from .sqe import StableHardness

for _value_type in (Profile, Transform):
    _value_type.__module__ = __name__

__all__ = ["PROFILES", "Profile", "StableHardness", "Transform"]

PROFILES: Mapping[str, Profile] = MappingProxyType({**eem.PROFILES, **sqe.PROFILES})
