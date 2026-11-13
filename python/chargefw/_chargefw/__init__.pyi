"""Private native extension declarations for ChargeFW.

Keep this stub package aligned with the registrations in ``python/src/``.
"""

from . import calculation as calculation
from . import core as core
from . import io as io
from . import methods as methods
from . import parameterization as parameterization
from . import parameters as parameters

def version() -> str: ...
