"""Presenting classes under the public module that exports them."""

from typing import get_type_hints


def set_public_module(module: str, *types: type) -> None:
    """Report `module` as the module of each type, keeping its annotations resolvable.

    typing.get_type_hints() resolves string annotations in the module named by __module__, so they
    are resolved in the defining module before it is replaced.
    """

    for value_type in types:
        value_type.__annotations__ = get_type_hints(value_type)
        value_type.__module__ = module
