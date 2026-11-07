"""Parameter templates, evaluation of candidate parameter values, and fitting."""

import json
from pathlib import Path
from unittest.mock import patch

import chargefw
import numpy as np
import pytest
from chargefw.parameterisation import (
    EvaluationError,
    Evaluator,
    Slot,
    fit,
    make_parameter_template,
)
from chargefw.parameterisation.profiles import PROFILES, StableHardness
from numpy.typing import ArrayLike, NDArray


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


def test_native_normal_equations_match_the_jacobian() -> None:
    evaluator = eem_evaluator(water())
    native = evaluator._native
    assert native._exact_derivatives()
    references = np.array([-0.5, 0.3, 0.2])
    weights = np.array([1.0, 2.0, 0.5])

    payload = native._normal_equations(EEM_WATER_VALUES, references, weights)
    assert payload["status"] == "success"
    charges = evaluator.charges(EEM_WATER_VALUES)
    np.testing.assert_array_equal(payload["charges"], charges)
    residuals = weights * (charges - references)
    assert payload["sum_of_squares"] == pytest.approx(residuals @ residuals)
    weighted = weights[:, np.newaxis] * evaluator.jacobian(EEM_WATER_VALUES)
    gauss_newton = np.asarray(payload["gauss_newton"]).reshape(5, 5)
    np.testing.assert_allclose(gauss_newton, weighted.T @ weighted, atol=1e-5)
    np.testing.assert_allclose(np.asarray(payload["gradient"]), weighted.T @ residuals, atol=1e-5)

    # A failing calculation is a status; malformed arguments raise.
    failed = native._normal_equations(np.array([0.0, 2.0, 1.0, 3.0, 1.0]), references, weights)
    assert failed["status"] == "numerical_failure" and failed["charges"] is None
    assert failed["failure_message"] is not None and "kappa" in failed["failure_message"]
    with pytest.raises(ValueError):
        native._normal_equations(EEM_WATER_VALUES, references[:2], weights)


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


def training_evaluator() -> Evaluator:
    molecules = chargefw.MoleculeCollection((water(), formaldehyde()))
    template = make_parameter_template(molecules, "eem", atom_types="plain")
    plan = chargefw.assess(
        molecules, method="eem", parameter_set=template, execution="full"
    ).default_plan
    assert plan is not None
    return Evaluator(plan, threads=1)


# EEM values for kappa, then hydrogen, carbon, and oxygen (A, B), from which references are made.
TRUE_VALUES = np.array([0.5, 2.0, 1.0, 2.5, 1.2, 3.0, 1.1])
FIT_START = np.array([0.5, 2.0, 0.9, 2.2, 1.0, 2.6, 1.0])
# kappa and hydrogen A stay fixed: EEM charges do not change when every A shifts by the same amount
# or when A, B, and kappa scale together.
FIT_FREE = [2, 3, 4, 5, 6]


def test_fit_reproduces_references_made_from_known_values() -> None:
    evaluator = training_evaluator()
    references = evaluator.charges(TRUE_VALUES)

    result = fit(evaluator, references, start=FIT_START, free=FIT_FREE)

    assert result.success, result.message
    assert result.rmsd < 1e-8
    assert result.free == tuple(FIT_FREE)
    assert result.values[0] == FIT_START[0] and result.values[1] == FIT_START[1]
    np.testing.assert_array_equal(result.start, FIT_START)
    assert result.failed_evaluations == 0
    np.testing.assert_allclose(evaluator.charges(result.values), references, atol=1e-8)


def test_fit_takes_start_by_key_and_references_per_assignment() -> None:
    evaluator = training_evaluator()
    per_assignment = [evaluator.charges(TRUE_VALUES)[:3], evaluator.charges(TRUE_VALUES)[3:]]
    start = {(slot.key, slot.name): value for slot, value in zip(evaluator.layout, FIT_START)}

    result = fit(evaluator, per_assignment, start=start, free=FIT_FREE)
    assert result.rmsd < 1e-8

    oxygen = chargefw.AtomParameterKey(element="O", classifier="plain", type="*")
    with pytest.raises(KeyError):
        fit(evaluator, per_assignment, start={(oxygen, "C"): 1.0}, free=FIT_FREE)
    with pytest.raises(ValueError, match="one per calculated charge"):
        fit(evaluator, per_assignment[0], start=FIT_START, free=FIT_FREE)
    # The right total, split at the wrong place, would pair references with the wrong atoms.
    flat = evaluator.charges(TRUE_VALUES)
    with pytest.raises(ValueError, match="assignment 0 must have 3 values"):
        fit(evaluator, [flat[:4], flat[4:]], start=FIT_START, free=FIT_FREE)
    with pytest.raises(ValueError, match="2 assignments"):
        fit(evaluator, [flat], start=FIT_START, free=FIT_FREE)
    with pytest.raises(ValueError, match="weights for assignment 1"):
        fit(evaluator, flat, weights=[np.ones(3), np.ones(3)], start=FIT_START, free=FIT_FREE)


def test_fit_checks_the_start_before_fitting() -> None:
    evaluator = training_evaluator()
    references = evaluator.charges(TRUE_VALUES)
    lower = np.full(7, -np.inf)
    upper = np.full(7, np.inf)

    upper[5] = 2.0
    with pytest.raises(ValueError, match="slot 5 \\(atom A of"):
        fit(evaluator, references, start=FIT_START, free=FIT_FREE, bounds=(lower, upper))
    upper[5] = FIT_START[5]
    with pytest.raises(ValueError, match="strictly inside"):
        fit(evaluator, references, start=FIT_START, free=FIT_FREE, bounds=(lower, upper))
    # So close to a bound that SciPy would move the start itself.
    upper[5] = FIT_START[5] + 1e-12
    with pytest.raises(ValueError, match="strictly inside"):
        fit(evaluator, references, start=FIT_START, free=FIT_FREE, bounds=(lower, upper))

    zero_kappa = FIT_START.copy()
    zero_kappa[0] = 0.0
    with pytest.raises(EvaluationError, match="at the start"):
        fit(evaluator, references, start=zero_kappa, free=FIT_FREE)

    # Equal electronegativities give zero charges for neutral molecules, which no hardness changes.
    equal = FIT_START.copy()
    equal[[1, 3, 5]] = 2.0
    with pytest.raises(ValueError, match="slot 2 \\(atom B of"):
        fit(evaluator, references, start=equal, free=FIT_FREE)

    with pytest.raises(ValueError, match="at least one free slot"):
        fit(evaluator, references, start=FIT_START, free=[])
    with pytest.raises(ValueError, match="must not repeat"):
        fit(evaluator, references, start=FIT_START, free=[5, 5])


class FailingEvaluator(Evaluator):
    """Fails every calculation with an oxygen electronegativity above 2.9."""

    def charges(self, values: ArrayLike) -> NDArray[np.float64]:
        vector = np.asarray(values, dtype=np.float64)
        if vector[5] > 2.9:
            raise EvaluationError("oxygen electronegativity above 2.9")
        return super().charges(vector)


def test_fit_rejects_trial_values_whose_calculation_fails() -> None:
    molecules = chargefw.MoleculeCollection((water(), formaldehyde()))
    template = make_parameter_template(molecules, "eem", atom_types="plain")
    plan = chargefw.assess(
        molecules, method="eem", parameter_set=template, execution="full"
    ).default_plan
    assert plan is not None
    evaluator = FailingEvaluator(plan, threads=1)
    references = training_evaluator().charges(TRUE_VALUES)

    result = fit(evaluator, references, start=FIT_START, free=FIT_FREE)

    # The true oxygen value, 3.0, cannot be calculated, so the fit settles below 2.9.
    assert result.failed_evaluations > 0
    assert result.values[5] <= 2.9


def test_evaluator_jacobian_steps_the_other_way_when_a_step_fails() -> None:
    molecules = chargefw.MoleculeCollection((water(), formaldehyde()))
    template = make_parameter_template(molecules, "eem", atom_types="plain")
    plan = chargefw.assess(
        molecules, method="eem", parameter_set=template, execution="full"
    ).default_plan
    assert plan is not None
    near_failure = FIT_START.copy()
    near_failure[5] = 2.9

    failing = FailingEvaluator(plan, threads=1).jacobian(near_failure, slots=[5])
    expected = training_evaluator().jacobian(near_failure, slots=[5])
    np.testing.assert_allclose(failing, expected, rtol=1e-4)


def test_fitted_values_become_a_parameter_set_to_calculate_and_write(tmp_path: Path) -> None:
    molecules = chargefw.MoleculeCollection((water(), formaldehyde()))
    evaluator = training_evaluator()
    references = evaluator.charges(TRUE_VALUES)
    result = fit(evaluator, references, start=FIT_START, free=FIT_FREE)

    fitted = result.parameter_set(id="water-fit", name="Water fit", notes="Fitted in a test")
    assert (fitted.id, fitted.method, fitted.name, fitted.notes) == (
        "water-fit",
        "eem",
        "Water fit",
        "Fitted in a test",
    )
    calculated = chargefw.calculate(molecules, method="eem", parameter_set=fitted)
    charges = np.concatenate([assignment.values for assignment in calculated.assignments])
    np.testing.assert_allclose(charges, references, atol=1e-8)

    # Written and read back, the set gives exactly the same charges.
    path = tmp_path / "water-fit.json"
    chargefw.write_parameter_set(fitted, path)
    reread = chargefw.read_parameter_set(path)
    assert (reread.id, reread.name, reread.atom_keys) == (fitted.id, fitted.name, fitted.atom_keys)
    recalculated = chargefw.calculate(molecules, method="eem", parameter_set=reread)
    np.testing.assert_array_equal(
        np.concatenate([assignment.values for assignment in recalculated.assignments]), charges
    )

    with pytest.raises(ValueError, match="bundled"):
        result.parameter_set(id="EEM_Baek1991")


def test_write_parameter_set_writes_bundled_sets(tmp_path: Path) -> None:
    path = tmp_path / "baek.json"
    chargefw.write_parameter_set(chargefw.parameter_sets["EEM_Baek1991"], path)

    document = json.loads(path.read_text(encoding="utf-8"))
    assert document["metadata"]["id"] == "EEM_Baek1991"
    assert document["metadata"]["method"] == "eem"
    assert document["common"]["names"] == ["kappa"]


def template_evaluator(method: str, atom_types: str, bond_types: str = "plain") -> Evaluator:
    molecules = chargefw.MoleculeCollection((water(), formaldehyde()))
    template = make_parameter_template(
        molecules,
        method,
        atom_types=atom_types,  # type: ignore[arg-type]
        bond_types=bond_types,  # type: ignore[arg-type]
    )
    plan = chargefw.assess(
        molecules, method=method, parameter_set=template, execution="full"
    ).default_plan
    assert plan is not None
    return Evaluator(plan, threads=1)


def test_eem_profile_fits_a_template_from_the_electronegativities_alone() -> None:
    evaluator = template_evaluator("eem", "plain")
    references = training_evaluator().charges(TRUE_VALUES)
    carbon = chargefw.AtomParameterKey(element="C", classifier="plain", type="*")
    oxygen = chargefw.AtomParameterKey(element="O", classifier="plain", type="*")

    # Hydrogen A is the fixed slot, so it needs no value; kappa and B start from the profile.
    result = fit(evaluator, references, start={(carbon, "A"): 0.5, (oxygen, "A"): 1.5})

    assert result.rmsd < 1e-8, result.message
    assert result.free == (2, 3, 4, 5, 6)
    assert result.values[0] == PROFILES["eem"].start["kappa"]
    np.testing.assert_array_equal(result.bounds[0][[2, 4, 6]], [1e-3, 1e-3, 1e-3])


def test_profile_reports_free_slots_without_a_start_value() -> None:
    evaluator = template_evaluator("eem", "plain")
    references = training_evaluator().charges(TRUE_VALUES)

    with pytest.raises(ValueError, match="no start value") as error:
        fit(evaluator, references)
    message = str(error.value)
    assert "slot 3 (atom A of" in message and "slot 5 (atom A of" in message
    assert "slot 1 " not in message

    # Without a profile, the template's zero values are the start, and EEM rejects kappa = 0.
    with pytest.raises(EvaluationError, match="at the start"):
        fit(evaluator, references, profile=None)


# Atom types H/1, C/2, O/1, O/2 (electronegativity, hardness, width), then bonds H-C, H-O, C-O.
SQE_TRUE_VALUES = np.array(
    [4.5, 13.0, 1.0, 5.5, 10.0, 1.0, 8.7, 13.0, 1.0, 8.0, 11.0, 1.0, 1.0, 1.5, 0.8]
)


def sqe_start(
    evaluator: Evaluator,
) -> dict[tuple[chargefw.AtomParameterKey | chargefw.BondParameterKey | None, str], float]:
    return {
        (slot.key, slot.name): value
        for slot, value in zip(evaluator.layout, SQE_TRUE_VALUES - 0.5, strict=True)
        if slot.name == "electronegativity"
    }


def test_sqe_profile_keeps_every_fitted_type_stable() -> None:
    evaluator = template_evaluator("sqe", "hbo", "bo")
    references = evaluator.charges(SQE_TRUE_VALUES)
    start = sqe_start(evaluator)

    result = fit(evaluator, references, start=start)

    assert result.rmsd < 1e-6, result.message
    assert 0 not in result.free
    # The stability condition, per atom type: sqrt(pi) * hardness * width >= 1, and kappa >= 0.01.
    hardness, width = result.values[[1, 4, 7, 10]], result.values[[2, 5, 8, 11]]
    assert np.all(np.sqrt(np.pi) * hardness * width >= 1.0 - 1e-12)
    assert np.all(result.values[12:] >= 0.01)

    soft = dict(start)
    soft[(evaluator.layout[1].key, "hardness")] = 0.5
    with pytest.raises(ValueError, match="slot 1 \\(atom hardness of .*, fitted as u"):
        fit(evaluator, references, start=soft)

    # Narrow widths with a large hardness are stable; a width then moves the charges through the
    # hardness that depends on it, so it is not insensitive.
    narrow = dict(start)
    for slot in evaluator.layout:
        if slot.name in ("hardness", "width"):
            narrow[(slot.key, slot.name)] = 13.0 if slot.name == "hardness" else 0.06
    assert fit(evaluator, references, start=narrow, max_evaluations=1).evaluations >= 1


def test_stable_hardness_pairs_each_entry_with_its_own_width(tmp_path: Path) -> None:
    # Two entries share a key; ChargeFW uses the first, and each has its own hardness and width.
    hydrogen = {"element": "H", "classifier": "plain", "type": "*"}
    oxygen = {"element": "O", "classifier": "plain", "type": "*"}
    document = {
        "metadata": {"id": "repeated-keys", "method": "sqe"},
        "atom": {
            "names": ["electronegativity", "hardness", "width"],
            "data": [
                {"key": hydrogen, "values": [4.5, 2.0, 0.1]},
                {"key": hydrogen, "values": [4.5, 2.0, 10.0]},
                {"key": oxygen, "values": [8.7, 13.0, 1.0]},
            ],
        },
        "bond": {
            "names": ["kappa"],
            "data": [
                {
                    "key": {
                        "atoms": [hydrogen, oxygen],
                        "bond": {"classifier": "plain", "type": "*"},
                    },
                    "values": [1.0],
                }
            ],
        },
    }
    path = tmp_path / "repeated.json"
    path.write_text(json.dumps(document))
    parameter_set = chargefw.read_parameter_set(path)
    plan = chargefw.assess(
        water(), method="sqe", parameter_set=parameter_set, execution="full"
    ).default_plan
    assert plan is not None
    evaluator = Evaluator(plan)
    names = [slot.name for slot in evaluator.layout]

    pairs = StableHardness().dependencies(evaluator.layout)

    assert pairs == {1: (2,), 4: (5,), 7: (8,)}, names
    # The first entry, with width 0.1, is outside the stable region: u < 0.
    assert StableHardness().to_internal(evaluator.values, evaluator.layout)[1] < 0.0


def test_sqe_profile_rejects_settings_its_transform_cannot_honour() -> None:
    evaluator = template_evaluator("sqe", "hbo", "bo")
    references = evaluator.charges(SQE_TRUE_VALUES)
    start = sqe_start(evaluator)

    lower = np.full(len(evaluator.layout), -np.inf)
    lower[1] = 1.0
    with pytest.raises(ValueError, match="unbounded.*slot 1 "):
        fit(evaluator, references, start=start, bounds=(lower, np.inf))

    # A fixed hardness with a free width would change the hardness anyway.
    free = [slot for slot in range(len(evaluator.layout)) if slot not in (0, 1)]
    with pytest.raises(ValueError, match="not free.*slot 1 "):
        fit(evaluator, references, start=start, free=free)


def test_fit_diagnostics_show_what_the_data_cannot_determine() -> None:
    evaluator = training_evaluator()
    references = evaluator.charges(TRUE_VALUES)

    # Carbon occurs once, in formaldehyde, so its A and B cannot be told apart.
    result = fit(evaluator, references, start=FIT_START, free=FIT_FREE)
    assert result.diagnostics.r2 == pytest.approx(1.0, abs=1e-12)
    assert result.diagnostics.rank == 4
    assert set(result.diagnostics.weakest) == {3, 4}
    assert result.diagnostics.singular_values.shape == (5,)

    # Freeing kappa and hydrogen A adds EEM's two invariances, which nothing can determine.
    everything = fit(evaluator, references, start=FIT_START, free=list(range(7)))
    assert everything.diagnostics.rank == 4


def test_fit_diagnostics_report_bounds_and_failed_calculations() -> None:
    evaluator = training_evaluator()
    references = evaluator.charges(TRUE_VALUES)
    lower = np.full(7, -np.inf)
    upper = np.full(7, np.inf)
    upper[6] = 1.05
    start = TRUE_VALUES.copy()
    start[6] = 1.0

    # The true oxygen B, 1.1, lies above the bound, and no other slot is free to compensate.
    bounded = fit(evaluator, references, start=start, free=[6], bounds=(lower, upper))
    assert bounded.diagnostics.at_bounds == (6,)
    assert bounded.values[6] == pytest.approx(1.05)

    molecules = chargefw.MoleculeCollection((water(), formaldehyde()))
    template = make_parameter_template(molecules, "eem", atom_types="plain")
    plan = chargefw.assess(
        molecules, method="eem", parameter_set=template, execution="full"
    ).default_plan
    assert plan is not None
    failing = fit(FailingEvaluator(plan, threads=1), references, start=FIT_START, free=FIT_FREE)
    assert failing.diagnostics.failed_messages == ("oxygen electronegativity above 2.9",)


def test_fit_diagnostics_cover_more_free_slots_than_charges() -> None:
    evaluator = eem_evaluator(water())
    references = evaluator.charges(EEM_WATER_VALUES)
    start = EEM_WATER_VALUES + np.array([0.0, 0.0, 0.1, 0.2, 0.1])

    # Three charges with a fixed total determine at most two combinations of the five slots.
    result = fit(evaluator, references, start=start, free=[0, 1, 2, 3, 4])

    assert result.diagnostics.singular_values.shape == (5,)
    assert result.diagnostics.rank == 2
    assert result.diagnostics.weakest


def test_fit_diagnostics_do_not_build_a_matrix_in_the_number_of_charges() -> None:
    evaluator = training_evaluator()
    references = evaluator.charges(TRUE_VALUES)

    # A full SVD of the tall Jacobian would build a square matrix in the number of charges.
    with patch.object(np.linalg, "svd", wraps=np.linalg.svd) as svd:
        fit(evaluator, references, start=FIT_START, free=FIT_FREE)

    assert svd.call_args.kwargs["full_matrices"] is False


def test_evaluator_counts_the_atoms_of_each_type() -> None:
    evaluator = training_evaluator()

    # kappa, then hydrogen (4 atoms), carbon (1), and oxygen (2), each with A and B.
    assert evaluator.matches == (None, 4, 4, 1, 1, 2, 2)


def test_fit_diagnostics_show_counts_and_errors_per_element() -> None:
    evaluator = training_evaluator()
    references = evaluator.charges(TRUE_VALUES)

    exact = fit(evaluator, references, start=FIT_START, free=FIT_FREE)
    # The weak direction lies on the one carbon atom, which the counts make visible.
    assert exact.diagnostics.matches == (4, 1, 1, 2, 2)
    assert set(exact.diagnostics.weakest) == {3, 4}
    by_element = exact.diagnostics.by_element
    assert list(by_element) == ["H", "C", "O"]
    assert [element.count for element in by_element.values()] == [4, 1, 2]
    assert all(element.rmsd < 1e-8 for element in by_element.values())
    # Four hydrogen charges vary; R² is defined. A single carbon charge leaves R² undefined.
    assert by_element["H"].r2 == pytest.approx(1.0, abs=1e-8)
    assert np.isnan(by_element["C"].r2)

    # Oxygen B held below its true value: the error is largest on oxygen.
    lower = np.full(7, -np.inf)
    upper = np.full(7, np.inf)
    upper[6] = 1.05
    start = TRUE_VALUES.copy()
    start[6] = 1.0
    bounded = fit(evaluator, references, start=start, free=[6], bounds=(lower, upper))
    errors = bounded.diagnostics.by_element
    assert max(errors, key=lambda element: errors[element].rmsd) == "O"
    assert errors["H"].r2 < 1.0 and errors["O"].r2 < 1.0


def test_matches_count_bonds_and_conformers_line_up_with_elements() -> None:
    # SQE counts bonds too: H-C twice, H-O twice, C-O once.
    sqe = template_evaluator("sqe", "hbo", "bo")
    bond_matches = [match for slot, match in zip(sqe.layout, sqe.matches) if slot.scope == "bond"]
    assert bond_matches == [2, 2, 1]

    # EEM uses coordinates, so water with two conformers gives two assignments of three charges.
    two_conformers = chargefw.Molecule(
        [8, 1, 1],
        bonds=np.array([[0, 1, 1], [0, 2, 1]], dtype=np.int32),
        coordinates=np.array(
            [
                [[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]],
                [[0.0, 0.0, 0.0], [0.98, 0.0, 0.0], [-0.26, 0.92, 0.0]],
            ]
        ),
    )
    evaluator = eem_evaluator(two_conformers)
    references = evaluator.charges(EEM_WATER_VALUES)
    result = fit(evaluator, references, start=EEM_WATER_VALUES + 0.1, free=[2, 3, 4])

    # Atoms are counted once per molecule, charges once per conformer.
    assert evaluator.matches == (None, 2, 2, 1, 1)
    by_element = result.diagnostics.by_element
    assert (by_element["H"].count, by_element["O"].count) == (4, 2)
