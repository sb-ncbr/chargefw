"""Focused smoke test for the installed package."""

import os
from pathlib import Path
from typing import get_type_hints

import chargefw
import chargefw.io.gemmi  # noqa: F401
import chargefw.io.rdkit  # noqa: F401
import chargefw.parameterization
import chargefw.parameterization.profiles


def test_version_and_typing_marker() -> None:
    assert isinstance(chargefw.__version__, str)
    assert chargefw.__version__
    expected_version = os.environ.get("CHARGEFW_EXPECTED_VERSION")
    if expected_version is not None:
        assert chargefw.__version__ == expected_version
    assert (Path(chargefw.__file__).parent / "py.typed").is_file()


def test_public_classes_report_a_public_module() -> None:
    # Renaming __module__ must keep annotations resolvable by typing.get_type_hints().
    for module in (
        chargefw,
        chargefw.calculation,
        chargefw.parameterization,
        chargefw.parameterization.profiles,
    ):
        for name in module.__all__:
            value = getattr(module, name)
            if isinstance(value, type):
                assert not any(part.startswith("_") for part in value.__module__.split(".")), (
                    module.__name__,
                    name,
                )
                get_type_hints(value)
