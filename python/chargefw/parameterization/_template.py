"""Parameter sets with an entry for every atom and bond type of some molecules."""

from __future__ import annotations

from collections.abc import Iterable

from .._chargefw import parameterization as _native_parameterization
from .._methods import Method
from .._parameters import ParameterSet, _parameter_set
from .._types import AtomClassifier, BondClassifier
from ..calculation import _as_collection, parameter_sets
from ..core import Molecule, MoleculeCollection


def make_parameter_template(
    molecules: Molecule | MoleculeCollection | Iterable[Molecule],
    method: str | Method,
    *,
    atom_types: AtomClassifier,
    bond_types: BondClassifier = "plain",
    id: str = "template",
) -> ParameterSet:
    """Make a parameter set with an entry for every atom and bond type of the molecules.

    Atom types are derived with the ``atom_types`` classifier and bond types with ``bond_types``,
    as described in the parameter-set format. Every entry has the parameter names the method needs,
    with all values zero. Every atom and bond of the molecules matches an entry. The ID must differ
    from the IDs of the bundled parameter sets.
    """

    if id in parameter_sets:
        raise ValueError(f"parameter-set ID {id!r} is the ID of a bundled parameter set")
    method_id = method.id if isinstance(method, Method) else method
    collection = _as_collection(molecules)
    native = _native_parameterization._make_parameter_template_catalog(
        collection._native_molecules, method_id, atom_types, bond_types, id
    )
    (descriptor,) = native._descriptors()
    return _parameter_set(descriptor, native)
