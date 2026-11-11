"""Focused smoke test for the installed package."""

import os
from pathlib import Path

import chargefw
import chargefw.io.gemmi  # noqa: F401
import chargefw.io.rdkit  # noqa: F401


def test_version_and_typing_marker() -> None:
    assert isinstance(chargefw.__version__, str)
    assert chargefw.__version__
    expected_version = os.environ.get("CHARGEFW_EXPECTED_VERSION")
    if expected_version is not None:
        assert chargefw.__version__ == expected_version
    assert (Path(chargefw.__file__).parent / "py.typed").is_file()


def test_public_classes_report_a_public_module() -> None:
    for name in chargefw.__all__:
        value = getattr(chargefw, name)
        if isinstance(value, type):
            assert not any(part.startswith("_") for part in value.__module__.split(".")), name
