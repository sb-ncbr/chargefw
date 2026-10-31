"""Conversion from upstream Gemmi objects to ChargeFW molecules."""

from __future__ import annotations

from collections.abc import Sequence
from importlib import import_module
from typing import Any, Literal, cast

from .._chargefw import adapters as _native_adapters
from ..calculation import CalculationResult
from ..core import MoleculeCollection
from . import (
    BondStrategy,
    ConformerSelection,
    RecordSelection,
    _native_results,
    dumps,
    parse,
)

try:
    import gemmi as _gemmi
except ModuleNotFoundError as error:
    if error.name != "gemmi":
        raise
    _gemmi = None  # type: ignore[assignment]


def _require_gemmi() -> Any:
    try:
        return import_module("gemmi")
    except ModuleNotFoundError as error:
        if error.name != "gemmi":
            raise
        raise ImportError(
            "chargefw.io.gemmi requires the optional Gemmi Python integration; "
            "install it with `pip install chargefw[gemmi]`"
        ) from error


# A gemmi.Structure keeps its links between residues (structure.connections) but neither the
# file's _chem_comp_bond rows nor the bond orders of the links, and writing it as mmCIF drops its
# PDB CONECT records. The regenerated document therefore has _struct_conn rows and no
# _chem_comp_bond, so "hybrid" yields template bonds plus the links; "explicit" would yield the
# links alone.
type StructureBondStrategy = Literal["none", "templates", "hybrid"]


def from_structure(
    structure: _gemmi.Structure,
    *,
    source_name: str = "",
    selection: RecordSelection = "all",
    bonds: StructureBondStrategy = "hybrid",
    conformers: ConformerSelection = "all",
) -> MoleculeCollection:
    """Convert an in-memory ``gemmi.Structure`` to a one-molecule collection.

    With ``bonds="hybrid"``, residues get their template bonds, and the covalent and disulfide links
    in ``structure.connections`` are added as single bonds. ``"templates"`` leaves the links out.
    Residues without a template, such as ligands and modified residues, get no bonds within them;
    the links still connect them to other residues. PDB ``CONECT`` records (``conect_map``) are not
    used. For the bonds the file records, use ``from_document()`` for mmCIF or
    ``chargefw.io.read()``.
    """

    gemmi = _require_gemmi()
    if not isinstance(structure, gemmi.Structure):
        raise TypeError("structure must be a gemmi.Structure")
    if bonds not in ("none", "templates", "hybrid"):
        raise ValueError(
            f"from_structure() does not support bonds={bonds!r}: it uses only residue templates "
            "and the links in structure.connections. For mmCIF, keep the document "
            "(gemmi.cif.read(path), or gemmi.read_structure(path, save_doc=document)) and use "
            "chargefw.io.gemmi.from_document(document); for PDB, use "
            f"chargefw.io.read(path, format='pdb', bonds={bonds!r})"
        )
    return parse(
        structure.make_mmcif_document().as_string(),
        format="mmcif",
        source_name=source_name,
        selection=selection,
        bonds=bonds,
        conformers=conformers,
    )


def from_document(
    document: _gemmi.cif.Document,
    *,
    source_name: str = "",
    selection: RecordSelection = "all",
    bonds: BondStrategy = "hybrid",
    conformers: ConformerSelection = "all",
) -> MoleculeCollection:
    """Convert a ``gemmi.cif.Document`` through ChargeFW's mmCIF parser."""

    gemmi = _require_gemmi()
    if not isinstance(document, gemmi.cif.Document):
        raise TypeError("document must be a gemmi.cif.Document")
    return parse(
        document.as_string(),
        format="mmcif",
        source_name=source_name,
        selection=selection,
        bonds=bonds,
        conformers=conformers,
    )


def to_document(
    result: CalculationResult | Sequence[CalculationResult],
) -> _gemmi.cif.Document:
    """Create a fresh Gemmi mmCIF document from one or more results for the same molecules."""

    gemmi = _require_gemmi()
    return cast("_gemmi.cif.Document", gemmi.cif.read_string(dumps(result, format="mmcif")))


def attach_charges(
    document: _gemmi.cif.Document,
    result: CalculationResult | Sequence[CalculationResult],
    *,
    overwrite: bool = False,
) -> None:
    """Attach charges to the unchanged Gemmi document used for the calculation input.

    Several results for the same molecules are attached as separate charge sets.
    """

    gemmi = _require_gemmi()
    if not isinstance(document, gemmi.cif.Document):
        raise TypeError("document must be a gemmi.cif.Document")
    results = _native_results(result)
    if not isinstance(overwrite, bool):
        raise TypeError("overwrite must be a bool")
    attached = gemmi.cif.read_string(
        _native_adapters._attach_mmcif(document.as_string(), results, overwrite)
    )
    document.clear()
    for block in attached:
        document.add_copied_block(block)


__all__ = ["from_structure", "from_document", "to_document", "attach_charges"]
