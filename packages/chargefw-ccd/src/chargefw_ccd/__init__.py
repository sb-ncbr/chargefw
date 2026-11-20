"""Atoms and bonds of the wwPDB Chemical Component Dictionary, as a data file for ChargeFW."""

from pathlib import Path

from ._version import __version__

__all__ = ["__version__", "path"]


def path() -> Path:
    """Returns the location of the data file."""
    return Path(__file__).parent / "components.txt.gz"
