"""Parameter templates and evaluation of candidate parameter values."""

import chargefw
import numpy as np
import pytest
from chargefw.parameterisation import EvaluationError, Evaluator, Slot, make_parameter_template


def water() -> chargefw.Molecule:
    return chargefw.Molecule(
        [8, 1, 1],
        bonds=np.array([[0, 1, 1], [0, 2, 1]], dtype=np.int32),
        coordinates=np.array([[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]]),
    )


def formaldehyde() -> chargefw.Molecule:
    return chargefw.Molecule(
        [6, 8, 1, 1],
        bonds=np.array([[0, 1, 2], [0, 2, 1], [3, 0, 1]], dtype=np.int32),
        coordinates=np.array(
            [[0.0, 0.0, 0.0], [1.21, 0.0, 0.0], [-0.55, 0.94, 0.0], [-0.55, -0.94, 0.0]]
        ),
    )


def eem_evaluator(molecule: chargefw.Molecule) -> Evaluator:
    template = make_parameter_template(molecule, "eem", atom_types="plain")
    plan = chargefw.assess(
        molecule, method="eem", parameter_set=template, execution="full"
    ).default_plan
    assert plan is not None
    return Evaluator(plan, threads=1)


# Slots of an EEM template for water: kappa, then hydrogen (A, B), then oxygen (A, B).
EEM_WATER_VALUES = np.array([0.5, 2.0, 1.0, 3.0, 1.0])


def test_parameter_template_has_the_types_of_its_molecules() -> None:
    molecules = chargefw.MoleculeCollection((water(), formaldehyde()))
    template = make_parameter_template(molecules, "sqe", atom_types="hbo", bond_types="bo")

    assert (template.id, template.method) == ("template", "sqe")
    assert [(key.element, key.type) for key in template.atom_keys] == [
        ("H", "1"),
        ("C", "2"),
        ("O", "1"),
        ("O", "2"),
    ]
    assert [
        (key.atoms[0].element, key.atoms[1].element, key.type) for key in template.bond_keys
    ] == [("H", "C", "1"), ("H", "O", "1"), ("C", "O", "2")]

    plan = chargefw.assess(
        molecules, method="sqe", parameter_set=template, execution="full"
    ).default_plan
    assert plan is not None
    evaluator = Evaluator(plan)
    # Three atom parameters per atom type and one bond parameter per bond type, all zero.
    assert len(evaluator.layout) == 4 * 3 + 3 * 1
    np.testing.assert_array_equal(evaluator.values, np.zeros(15))


def test_parameter_template_rejects_bundled_ids_and_methods_without_parameters() -> None:
    with pytest.raises(ValueError, match="bundled"):
        make_parameter_template(water(), "eem", atom_types="plain", id="EEM_Baek1991")
    with pytest.raises(ValueError):
        make_parameter_template(water(), "formal", atom_types="plain")
    with pytest.raises(ValueError):
        make_parameter_template(water(), "unknown", atom_types="plain")


def test_evaluator_layout_names_each_slot_by_key() -> None:
    evaluator = eem_evaluator(water())
    hydrogen = chargefw.AtomParameterKey(element="H", classifier="plain", type="*")
    oxygen = chargefw.AtomParameterKey(element="O", classifier="plain", type="*")

    assert evaluator.layout == (
        Slot("common", None, "kappa"),
        Slot("atom", hydrogen, "A"),
        Slot("atom", hydrogen, "B"),
        Slot("atom", oxygen, "A"),
        Slot("atom", oxygen, "B"),
    )


def test_evaluator_charges_follow_the_values() -> None:
    evaluator = eem_evaluator(water())

    charges = evaluator.charges(EEM_WATER_VALUES)
    assert charges.shape == (3,)
    assert charges[0] < 0.0 and charges.sum() == pytest.approx(0.0, abs=1e-12)

    equal = EEM_WATER_VALUES.copy()
    equal[3] = equal[1]
    np.testing.assert_allclose(evaluator.charges(equal), np.zeros(3), atol=1e-12)

    with pytest.raises(EvaluationError, match="kappa"):
        evaluator.charges(np.array([0.0, 2.0, 1.0, 3.0, 1.0]))
    with pytest.raises(ValueError):
        evaluator.charges(EEM_WATER_VALUES[:3])
    with pytest.raises(ValueError):
        evaluator.charges(np.array([0.5, np.nan, 1.0, 3.0, 1.0]))


def test_evaluator_jacobian_respects_eem_invariances_and_bounds() -> None:
    evaluator = eem_evaluator(water())
    jacobian = evaluator.jacobian(EEM_WATER_VALUES)
    assert jacobian.shape == (3, 5)

    # Adding the same constant to every A, or scaling A, B, and kappa together, changes no charge.
    shift = np.array([0.0, 1.0, 0.0, 1.0, 0.0])
    np.testing.assert_allclose(jacobian @ shift, np.zeros(3), atol=1e-6)
    np.testing.assert_allclose(jacobian @ EEM_WATER_VALUES, np.zeros(3), atol=1e-6)
    assert np.abs(jacobian[:, 3]).max() > 0.01

    selected = evaluator.jacobian(EEM_WATER_VALUES, slots=[3, 0])
    np.testing.assert_allclose(selected, jacobian[:, [3, 0]], atol=1e-6)

    # At an upper bound, the step goes down instead.
    upper = np.full(5, np.inf)
    upper[3] = EEM_WATER_VALUES[3]
    bounded = evaluator.jacobian(EEM_WATER_VALUES, slots=[3], bounds=(-upper, upper))
    np.testing.assert_allclose(bounded[:, 0], jacobian[:, 3], atol=1e-6)


def test_evaluator_jacobian_names_the_slot_whose_step_fails() -> None:
    evaluator = eem_evaluator(water())
    values = EEM_WATER_VALUES.copy()
    values[0] = 1e-12
    lower = np.full(5, -np.inf)
    upper = np.full(5, np.inf)
    upper[0] = values[0]

    # The bounds allow negative kappa, so the step goes down to a value EEM rejects.
    with pytest.raises(EvaluationError, match="slot 0 \\(common kappa\\)"):
        evaluator.jacobian(values, slots=[0], bounds=(lower, upper))


def test_evaluator_jacobian_steps_stay_within_narrow_bounds() -> None:
    evaluator = eem_evaluator(water())
    unbounded = evaluator.jacobian(EEM_WATER_VALUES, slots=[0])

    # An interval narrower than the step shortens it; the derivative stays the same.
    lower = np.full(5, -np.inf)
    upper = np.full(5, np.inf)
    lower[0] = EEM_WATER_VALUES[0] - 1e-10
    upper[0] = EEM_WATER_VALUES[0]
    narrow = evaluator.jacobian(EEM_WATER_VALUES, slots=[0], bounds=(lower, upper))
    np.testing.assert_allclose(narrow, unbounded, rtol=1e-4)

    # A tiny positive kappa in a tiny positive interval steps to its lower end, not below zero.
    values = EEM_WATER_VALUES.copy()
    values[0] = 1e-12
    lower[0], upper[0] = 5e-13, 1e-12
    assert np.all(np.isfinite(evaluator.jacobian(values, slots=[0], bounds=(lower, upper))))

    lower[0] = upper[0] = 1e-12
    with pytest.raises(ValueError, match="no room"):
        evaluator.jacobian(values, slots=[0], bounds=(lower, upper))
    with pytest.raises(ValueError, match="within the bounds"):
        evaluator.jacobian(EEM_WATER_VALUES, slots=[0], bounds=(lower, upper))

    # NaN is not "no limit": it would make every comparison with the bound false.
    lower[0], upper[0] = np.nan, np.inf
    with pytest.raises(ValueError, match="NaN"):
        evaluator.jacobian(EEM_WATER_VALUES, slots=[0], bounds=(lower, upper))


def test_evaluator_jacobian_takes_integer_slot_indices_only() -> None:
    evaluator = eem_evaluator(water())
    expected = evaluator.jacobian(EEM_WATER_VALUES, slots=[3])

    np.testing.assert_array_equal(
        evaluator.jacobian(EEM_WATER_VALUES, slots=[np.int64(3)]), expected
    )
    for invalid in (3.9, "3", True):
        with pytest.raises(TypeError):
            evaluator.jacobian(EEM_WATER_VALUES, slots=[invalid])  # type: ignore[list-item]
    for outside in (5, -1):
        with pytest.raises(IndexError):
            evaluator.jacobian(EEM_WATER_VALUES, slots=[outside])


def test_evaluator_needs_full_execution_and_parameters() -> None:
    molecule = water()
    template = make_parameter_template(molecule, "eem", atom_types="plain")
    cutoff = chargefw.assess(
        molecule, method="eem", parameter_set=template, execution="cutoff", radius=8.0
    ).default_plan
    assert cutoff is not None
    with pytest.raises(ValueError, match="full execution"):
        Evaluator(cutoff)

    formal = chargefw.assess(molecule, method="formal").default_plan
    assert formal is not None
    with pytest.raises(ValueError, match="no parameters"):
        Evaluator(formal)
