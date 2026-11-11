"""Python warnings for issues that can make calculated charges wrong.

Every issue stays recorded in its attribute, with its atoms. The warnings only make the issues that
can make charges wrong visible without looking: one warning per kind of issue and call, so that a
large batch does not flood the output.
"""

from collections.abc import Callable, Iterable
from typing import Protocol

from ._public import set_public_module


class ChargeFWWarning(UserWarning):
    """Base class of ChargeFW warnings."""


class InputWarning(ChargeFWWarning):
    """The input molecules have a problem that can make their charges wrong."""


class ChargeWarning(ChargeFWWarning):
    """Calculated charges are probably wrong."""


# Kinds that can make the charges wrong. Other kinds are only recorded.
IMPORT_WARNING_KINDS = frozenset({"unbonded_residue_atoms", "symmetry_connections_ignored"})
INPUT_WARNING_KINDS = frozenset({"hydrogens_missing", "atoms_too_close"})
CHARGE_WARNING_KINDS = frozenset({"charge_out_of_range"})


class _Issue(Protocol):
    @property
    def kind(self) -> str: ...

    @property
    def message(self) -> str: ...


def issue_warnings(
    issues: Iterable[tuple[int, _Issue]],
    *,
    kinds: frozenset[str],
    category: type[ChargeFWWarning],
    molecule_count: int,
    describe: Callable[[int], str],
    attribute: str,
) -> list[tuple[str, type[ChargeFWWarning]]]:
    """One warning text per kind for the issues of the given kinds, given as (molecule, issue).

    The public function the user called raises them, so that a warning turned into an error has a
    short traceback.
    """

    by_kind: dict[str, list[tuple[int, _Issue]]] = {}
    for molecule_index, issue in issues:
        if issue.kind in kinds:
            by_kind.setdefault(issue.kind, []).append((molecule_index, issue))
    result = []
    for kind, found in by_kind.items():
        first_index, first = found[0]
        example = f"{describe(first_index)}: {first.message.rstrip('.')}"
        affected = len({index for index, _ in found})
        if affected == 1:
            text = f"{example} ({kind}; see {attribute})"
        else:
            text = (
                f"{affected} of {molecule_count} molecules: {kind}, for example {example} "
                f"(see {attribute})"
            )
        result.append((text, category))
    return result


set_public_module("chargefw", ChargeFWWarning, InputWarning, ChargeWarning)
