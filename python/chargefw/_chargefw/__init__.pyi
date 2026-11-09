"""Private native extension declarations for ChargeFW.

Keep this stub package aligned with the registrations in ``python/src/``.
"""

from . import adapters as adapters
from . import calculation as calculation
from . import core as core
from . import methods as methods
from . import parameterisation as parameterisation
from . import parameters as parameters

def version() -> str: ...
