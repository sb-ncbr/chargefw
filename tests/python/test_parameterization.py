"""Parameter templates, evaluation of candidate parameter values, and fitting.

The tests check what the package promises, on molecules of a few atoms. They do not check how far an
optimization converges, nor parameter values along directions the data does not determine; fit
quality on real data is qualified outside the test suite.
"""

import json
from collections.abc import Callable
from pathlib import Path
from typing import get_type_hints
from unittest.mock import patch

import chargefw
import numpy as np
import pytest
from chargefw.parameterization import (
    ChargeDiagnostics,
    Diagnostics,
    EvaluationError,
    Evaluator,
    FitObserver,
    FitProgress,
    FitResult,
    Slot,
    fit,
    make_parameter_template,
)
from chargefw.parameterization._diagnostics import _spectrum_of_gram, _spectrum_of_jacobian
from chargefw.parameterization.profiles import PROFILES, Profile, StableHardness
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
# or when A, B, and kappa scale together. Carbon occurs once, so its B acts exactly as its A: the
# data cannot tell them apart.
FIT_FREE = [2, 3, 4, 5, 6]
# The slots these charges determine: carbon stays fixed too.
DETERMINED_FREE = [2, 5, 6]


def test_public_annotations_resolve() -> None:
    # Introspection such as typing.get_type_hints() finds every name the annotations use.
    for public in (
        Slot,
        Diagnostics,
        ChargeDiagnostics,
        FitProgress,
        FitResult,
        Evaluator.__init__,
        Evaluator.charges,
        fit,
        make_parameter_template,
    ):
        assert get_type_hints(public), public


def test_fit_recovers_the_charges_of_known_values() -> None:
    evaluator = training_evaluator()
    references = evaluator.charges(TRUE_VALUES)
    start = TRUE_VALUES.copy()
    start[DETERMINED_FREE] += 0.2

    result = fit(evaluator, references, start=start, free=DETERMINED_FREE)

    assert result.success, result.message
    assert result.free == tuple(DETERMINED_FREE)
    unchanged = [slot for slot in range(7) if slot not in DETERMINED_FREE]
    np.testing.assert_array_equal(result.values[unchanged], start[unchanged])
    np.testing.assert_array_equal(result.start, start)
    assert result.failed_evaluations == 0
    np.testing.assert_allclose(evaluator.charges(result.values), references, atol=1e-6)


def test_fit_takes_start_by_key_and_references_per_assignment() -> None:
    evaluator = training_evaluator()
    flat = evaluator.charges(TRUE_VALUES)
    per_assignment = [flat[:3], flat[3:]]
    start = {(slot.key, slot.name): value for slot, value in zip(evaluator.layout, FIT_START)}

    # The same fit, given by key and per assignment, takes the same steps.
    by_key = fit(evaluator, per_assignment, start=start, free=FIT_FREE, max_evaluations=5)
    by_vector = fit(evaluator, flat, start=FIT_START, free=FIT_FREE, max_evaluations=5)
    np.testing.assert_array_equal(by_key.start, FIT_START)
    np.testing.assert_array_equal(by_key.values, by_vector.values)

    oxygen = chargefw.AtomParameterKey(element="O", classifier="plain", type="*")
    with pytest.raises(KeyError):
        fit(evaluator, per_assignment, start={(oxygen, "C"): 1.0}, free=FIT_FREE)
    with pytest.raises(ValueError, match="one per calculated charge"):
        fit(evaluator, per_assignment[0], start=FIT_START, free=FIT_FREE)
    # The right total, split at the wrong place, would pair references with the wrong atoms.
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
    # A start on a bound, as an earlier fit can leave it, or beyond it by rounding, is accepted with
    # exact derivatives and with finite differences.
    differences = training_evaluator()
    differences._exact = False
    for path in (evaluator, differences):
        for bound in (FIT_START[5], FIT_START[5] - 1e-12):
            upper[5] = bound
            result = fit(
                path,
                references,
                start=FIT_START,
                free=FIT_FREE,
                bounds=(lower, upper),
                max_evaluations=3,
            )
            assert result.start[5] == FIT_START[5] and result.values[5] <= bound
    upper[5] = FIT_START[5] - 1e-6
    with pytest.raises(ValueError, match="within the bounds"):
        fit(evaluator, references, start=FIT_START, free=FIT_FREE, bounds=(lower, upper))
    upper[5] = np.inf

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

    def _normal_equations(
        self,
        values: NDArray[np.float64],
        references: NDArray[np.float64],
        weights: NDArray[np.float64],
    ) -> tuple[NDArray[np.float64], float, NDArray[np.float64], NDArray[np.float64]]:
        if values[5] > 2.9:
            raise EvaluationError("oxygen electronegativity above 2.9")
        return super()._normal_equations(values, references, weights)


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
    result = fit(evaluator, references, start=FIT_START, free=FIT_FREE, max_evaluations=5)

    fitted = result.parameter_set(id="water-fit", name="Water fit", notes="Fitted in a test")
    assert (fitted.id, fitted.method, fitted.name, fitted.notes) == (
        "water-fit",
        "eem",
        "Water fit",
        "Fitted in a test",
    )
    calculated = chargefw.calculate(molecules, method="eem", parameter_set=fitted)
    charges = np.concatenate([assignment.values for assignment in calculated.assignments])
    np.testing.assert_allclose(charges, evaluator.charges(result.values), atol=1e-12)

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


def test_eem_profile_fits_a_template_without_a_start() -> None:
    evaluator = template_evaluator("eem", "plain")
    references = training_evaluator().charges(TRUE_VALUES)

    # Hydrogen A is the fixed slot; kappa and B start from the profile, and the free A from their
    # best values for that start.
    result = fit(evaluator, references, max_evaluations=5)

    assert result.free == (2, 3, 4, 5, 6)
    profile = PROFILES["eem"].start
    np.testing.assert_array_equal(
        result.start[[0, 1, 2, 4, 6]],
        [profile["kappa"], profile["A"], profile["B"], profile["B"], profile["B"]],
    )
    assert result.values[0] == profile["kappa"] and result.values[1] == profile["A"]
    np.testing.assert_array_equal(result.bounds[0][[2, 4, 6]], [1e-3, 1e-3, 1e-3])


def test_profile_starts_linear_parameters_from_their_best_values() -> None:
    evaluator = template_evaluator("sqeqp", "hbo", "bo")
    references = training_evaluator().charges(TRUE_VALUES)
    weights = np.ones(references.size)
    linear = [
        slot
        for slot, entry in enumerate(evaluator.layout)
        if entry.name in ("electronegativity", "q0") and entry.key != evaluator.layout[0].key
    ]

    result = fit(evaluator, references, max_evaluations=1)

    # At the reported start, no change of the free linear slots lowers the sum of squares.
    assert set(linear) <= set(result.free)
    _, total, gram, gradient = evaluator._normal_equations(result.start, references, weights)
    scale = np.sqrt(np.diag(gram)[linear] * total)
    assert np.all(np.abs(gradient[linear]) <= 1e-8 * scale)
    # So the charges at the start do not depend on a start given for them, with exact derivatives
    # and with finite differences; values the data does not determine keep their start.
    charges = evaluator.charges(result.start)
    given = {(evaluator.layout[slot].key, evaluator.layout[slot].name): 3.0 for slot in linear}
    from_given = fit(evaluator, references, start=given, max_evaluations=1).start
    np.testing.assert_allclose(evaluator.charges(from_given), charges, atol=1e-10)
    differences = template_evaluator("sqeqp", "hbo", "bo")
    differences._exact = False
    from_differences = fit(differences, references, max_evaluations=1).start
    np.testing.assert_allclose(evaluator.charges(from_differences), charges, atol=1e-6)

    # A bounded linear slot keeps its start, and the others still fit around it.
    lower = np.full(len(evaluator.layout), -np.inf)
    upper = np.full(len(evaluator.layout), np.inf)
    lower[linear[0]], upper[linear[0]] = -10.0, 10.0
    bounded = fit(evaluator, references, start=given, bounds=(lower, upper), max_evaluations=1)
    assert bounded.start[linear[0]] == 3.0
    assert np.all(bounded.start[linear[1:]] != 3.0)


def test_profile_reports_free_slots_without_a_start_value() -> None:
    evaluator = template_evaluator("eem", "plain")
    references = training_evaluator().charges(TRUE_VALUES)
    without_a = Profile(start={"kappa": 1.0, "B": 2.0}, fixed=("kappa", "A"))

    with pytest.raises(ValueError, match="no start value") as error:
        fit(evaluator, references, profile=without_a)
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

    result = fit(evaluator, references, start=start, max_evaluations=20)

    assert 0 not in result.free
    # The stability condition, per atom type: sqrt(pi) * hardness * width >= 1, and kappa >= 0.01.
    hardness, width = result.values[[1, 4, 7, 10]], result.values[[2, 5, 8, 11]]
    assert np.all(np.sqrt(np.pi) * hardness * width >= 1.0 - 1e-12)
    assert np.all(result.values[12:] >= 0.01)

    # Carbon (C/2) is not terminal, so its hardness is fitted and must start stable.
    soft = dict(start)
    soft[(evaluator.layout[4].key, "hardness")] = 0.5
    with pytest.raises(ValueError, match="slot 4 \\(atom hardness of .*, fitted as u"):
        fit(evaluator, references, start=soft)

    # Narrow widths with a large hardness are stable; a width then moves the charges through the
    # hardness that depends on it, so it is not insensitive.
    narrow = dict(start)
    for slot in evaluator.layout:
        if slot.name in ("hardness", "width"):
            narrow[(slot.key, slot.name)] = 13.0 if slot.name == "hardness" else 0.06
    assert fit(evaluator, references, start=narrow, max_evaluations=1).evaluations >= 1


def test_exact_and_forward_difference_fits_agree() -> None:
    evaluator = training_evaluator()
    differences = training_evaluator()
    differences._exact = False
    references = evaluator.charges(TRUE_VALUES)
    start = TRUE_VALUES.copy()
    start[DETERMINED_FREE] += 0.2

    # Both paths recover the charges from the slots the data determines.
    for path in (evaluator, differences):
        result = fit(path, references, start=start, free=DETERMINED_FREE)
        np.testing.assert_allclose(evaluator.charges(result.values), references, atol=1e-6)

    # At the same values, the diagnostics from the normal equations match the SVD of the Jacobian.
    weights = np.linspace(0.5, 1.5, references.size)
    _, _, gram, _ = evaluator._normal_equations(start, references, weights)
    jacobian = weights[:, np.newaxis] * differences.jacobian(start)
    from_gram, _ = _spectrum_of_gram(gram)
    from_jacobian, _ = _spectrum_of_jacobian(jacobian)
    resolved = from_jacobian > 1e-3 * from_jacobian[0]
    np.testing.assert_allclose(from_gram[resolved], from_jacobian[resolved], rtol=1e-4)


def test_sqe_profile_holds_the_hardness_of_terminal_types_at_its_limit() -> None:
    # In water and formaldehyde, every H/1 and O/2 atom has one bond: their hardness acts only
    # together with the kappa of their bonds. The profile holds it at the stability limit, u = 0,
    # and still fits the width it depends on.
    evaluator = template_evaluator("sqe", "hbo", "bo")
    differences = template_evaluator("sqe", "hbo", "bo")
    differences._exact = False
    references = evaluator.charges(SQE_TRUE_VALUES)
    terminal_hardness = [1, 10]  # H/1 and O/2

    for path in (evaluator, differences):
        result = fit(path, references, start=sqe_start(path), max_evaluations=10)

        assert not set(terminal_hardness) & set(result.free)
        assert {2, 11} <= set(result.free)  # their widths
        for values in (result.start, result.values):
            hardness, width = values[terminal_hardness], values[[2, 11]]
            np.testing.assert_allclose(np.sqrt(np.pi) * hardness * width, 1.0, rtol=1e-12)

    # A start that already fits keeps its charges: the kappa of the terminal types' bonds take up
    # what their hardness gives away.
    for path in (evaluator, differences):
        warm = fit(path, references, start=SQE_TRUE_VALUES, max_evaluations=1)
        np.testing.assert_allclose(path.charges(warm.start), references, atol=1e-10)
        assert warm.rmsd < 1e-8

    # Fitting those slots yourself is still possible.
    every = fit(
        evaluator,
        references,
        start=sqe_start(evaluator),
        free=range(1, len(evaluator.layout)),
        max_evaluations=3,
    )
    assert set(terminal_hardness) <= set(every.free)


def test_sqe_profile_holds_no_hardness_whose_kappa_other_types_share(tmp_path: Path) -> None:
    # One kappa entry for every bond: the C-H and O-H bonds end in a terminal hydrogen, the C=O bond
    # does not, so no kappa change can take up a change of the hydrogen hardness alone.
    atom = {element: {"element": element, "classifier": "plain", "type": "*"} for element in "HCO"}
    anything = {"element": "*", "classifier": "plain", "type": "*"}
    document = {
        "metadata": {"id": "shared-kappa", "method": "sqe"},
        "atom": {
            "names": ["electronegativity", "hardness", "width"],
            "data": [
                {"key": atom["H"], "values": [4.5, 13.0, 1.0]},
                {"key": atom["C"], "values": [5.5, 10.0, 1.0]},
                {"key": atom["O"], "values": [8.7, 13.0, 1.0]},
            ],
        },
        "bond": {
            "names": ["kappa"],
            "data": [
                {
                    "key": {
                        "atoms": [anything, anything],
                        "bond": {"classifier": "plain", "type": "*"},
                    },
                    "values": [1.0],
                }
            ],
        },
    }
    path = tmp_path / "shared.json"
    path.write_text(json.dumps(document))
    molecules = chargefw.MoleculeCollection((water(), formaldehyde()))
    plan = chargefw.assess(
        molecules,
        method="sqe",
        parameter_set=chargefw.read_parameter_set(path),
        execution="full",
    ).default_plan
    assert plan is not None
    evaluator = Evaluator(plan, threads=1)

    result = fit(evaluator, evaluator.charges(evaluator.values), max_evaluations=1)

    hardness = [index for index, slot in enumerate(evaluator.layout) if slot.name == "hardness"]
    assert set(hardness) <= set(result.free)


def test_sqe_profile_holds_no_hardness_of_an_atom_bonded_to_a_terminal_atom() -> None:
    # In HCl, both atoms have one bond, so no single kappa change takes up a change of the hydrogen
    # hardness there, although water's hydrogens alone would allow it.
    hydrogen_chloride = chargefw.Molecule(
        [17, 1],
        bonds=np.array([[0, 1, 1]], dtype=np.int32),
        coordinates=np.array([[0.0, 0.0, 0.0], [1.27, 0.0, 0.0]]),
    )
    molecules = chargefw.MoleculeCollection((water(), hydrogen_chloride))
    template = make_parameter_template(molecules, "sqe", atom_types="plain")
    plan = chargefw.assess(
        molecules, method="sqe", parameter_set=template, execution="full"
    ).default_plan
    assert plan is not None
    evaluator = Evaluator(plan, threads=1)
    pauling = {"H": 2.20, "O": 3.44, "Cl": 3.16}
    start: dict[tuple[chargefw.AtomParameterKey | chargefw.BondParameterKey | None, str], float] = {
        (slot.key, slot.name): pauling[slot.key.element]
        for slot in evaluator.layout
        if slot.name == "electronegativity" and isinstance(slot.key, chargefw.AtomParameterKey)
    }
    references = np.array([-0.6, 0.3, 0.3, -0.2, 0.2])

    result = fit(evaluator, references, start=start, max_evaluations=1)

    hardness = [index for index, slot in enumerate(evaluator.layout) if slot.name == "hardness"]
    assert set(hardness) <= set(result.free)


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


def test_sqe_profile_rejects_settings_its_transform_cannot_honor() -> None:
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
    result = fit(evaluator, references, start=FIT_START, free=FIT_FREE, max_evaluations=5)
    assert result.diagnostics.rank == 4
    assert set(result.diagnostics.weakest) == {3, 4}
    assert result.diagnostics.singular_values.shape == (5,)

    # Freeing kappa and hydrogen A adds EEM's two invariances, which nothing can determine.
    everything = fit(evaluator, references, start=FIT_START, free=list(range(7)), max_evaluations=5)
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

    # With exact derivatives, the diagnostics come from the normal equations, without a Jacobian.
    with patch.object(np.linalg, "svd", wraps=np.linalg.svd) as svd:
        fit(evaluator, references, start=FIT_START, free=FIT_FREE, max_evaluations=5)
    assert not svd.called

    # With forward differences, a full SVD of the tall Jacobian would build a square matrix in the
    # number of charges.
    evaluator._exact = False
    with patch.object(np.linalg, "svd", wraps=np.linalg.svd) as svd:
        fit(evaluator, references, start=FIT_START, free=FIT_FREE, max_evaluations=5)
    assert svd.call_args.kwargs["full_matrices"] is False


def test_evaluator_counts_the_atoms_of_each_type() -> None:
    evaluator = training_evaluator()

    # kappa, then hydrogen (4 atoms), carbon (1), and oxygen (2), each with A and B.
    assert evaluator.matches == (None, 4, 4, 1, 1, 2, 2)


def test_fit_diagnostics_show_counts_and_errors_per_element() -> None:
    evaluator = training_evaluator()
    references = evaluator.charges(TRUE_VALUES)

    result = fit(evaluator, references, start=FIT_START, free=FIT_FREE, max_evaluations=5)
    # The undetermined direction lies on the one carbon atom, which the counts make visible.
    assert result.diagnostics.matches == (4, 1, 1, 2, 2)
    assert set(result.diagnostics.weakest) == {3, 4}
    by_element = result.diagnostics.by_element
    assert list(by_element) == ["H", "C", "O"]
    assert [element.count for element in by_element.values()] == [4, 1, 2]
    # The errors are those of the charges at the fitted values.
    charges = evaluator.charges(result.values)
    assert result.rmsd == pytest.approx(np.sqrt(np.mean((charges - references) ** 2)))
    assert result.diagnostics.r2 == pytest.approx(
        1.0 - np.sum((charges - references) ** 2) / np.sum((references - references.mean()) ** 2)
    )
    hydrogen = np.array([1, 2, 5, 6])
    assert by_element["H"].rmsd == pytest.approx(
        np.sqrt(np.mean((charges[hydrogen] - references[hydrogen]) ** 2))
    )
    assert by_element["H"].r2 == pytest.approx(
        1.0
        - np.sum((charges[hydrogen] - references[hydrogen]) ** 2)
        / np.sum((references[hydrogen] - references[hydrogen].mean()) ** 2)
    )
    # A single carbon charge leaves R² undefined.
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
    by_type = {key.element: value for key, value in result.diagnostics.by_atom_type.items()}
    assert (by_type["H"].count, by_type["O"].count) == (4, 2)


def test_fit_diagnostics_separate_atom_types_of_one_element() -> None:
    evaluator = template_evaluator("sqe", "hbo", "bo")
    references = evaluator.charges(SQE_TRUE_VALUES) + 0.01

    result = fit(evaluator, references, start=SQE_TRUE_VALUES, max_evaluations=3)

    # Water oxygen has single bonds only, formaldehyde oxygen a double bond: one element, two types.
    by_type = result.diagnostics.by_atom_type
    oxygen = {key.type: diagnostics for key, diagnostics in by_type.items() if key.element == "O"}
    assert set(oxygen) == {"1", "2"}
    assert (oxygen["1"].count, oxygen["2"].count) == (1, 1)
    assert result.diagnostics.by_element["O"].count == 2
    charges = evaluator.charges(result.values)
    # Water's oxygen comes first among the charges, formaldehyde's fifth.
    assert oxygen["1"].mean_error == pytest.approx(charges[0] - references[0])
    assert oxygen["2"].mean_error == pytest.approx(charges[4] - references[4])


class RecordingFitObserver(FitObserver):
    def __init__(self, stop_when: Callable[[FitProgress], bool] = lambda _: False) -> None:
        self.events: list[FitProgress] = []
        self._stop_when = stop_when

    def on_evaluation(self, progress: FitProgress) -> None:
        self.events.append(progress)

    def cancelled(self) -> bool:
        return self._stop_when(self.events[-1])


def test_observer_follows_every_trial_point_and_the_best_so_far() -> None:
    evaluator = training_evaluator()
    references = evaluator.charges(TRUE_VALUES)
    weights = np.array([1.0, 2.0, 2.0, 1.0, 3.0, 0.5, 0.5])
    observer = RecordingFitObserver()

    result = fit(
        evaluator,
        references,
        weights=weights,
        start=FIT_START,
        free=FIT_FREE,
        max_evaluations=8,
        observer=observer,
    )

    events = observer.events
    assert [event.evaluation for event in events] == list(range(1, len(events) + 1))
    assert len(events) == result.evaluations
    # The first trial point is the start.
    start_residuals = weights * (evaluator.charges(FIT_START) - references)
    assert events[0].objective == pytest.approx(start_residuals @ start_residuals)
    # The best point is the one with the lowest weighted objective, which the RMSD need not follow.
    for index, event in enumerate(events):
        lowest = min(range(index + 1), key=lambda earlier: events[earlier].objective)
        assert event.best_objective == events[lowest].objective
        assert event.best_rmsd == events[lowest].rmsd
    np.testing.assert_allclose(events[-1].best_values, result.values)
    assert all(
        later.elapsed_seconds >= earlier.elapsed_seconds
        for earlier, later in zip(events, events[1:])
    )


def test_cancelled_fit_returns_the_best_values_as_a_checkpoint() -> None:
    evaluator = training_evaluator()
    references = evaluator.charges(TRUE_VALUES)
    observer = RecordingFitObserver(lambda progress: progress.evaluation == 3)

    result = fit(evaluator, references, start=FIT_START, free=FIT_FREE, observer=observer)

    assert not result.success and "cancelled" in result.message
    assert result.evaluations == 3
    np.testing.assert_array_equal(result.values, observer.events[-1].best_values)
    # The best values so far, written as a parameter set, start another fit where this one stopped.
    checkpoint = observer.events[-1].parameter_set(id="checkpoint")
    resumed = fit(evaluator, references, start=checkpoint, free=FIT_FREE, max_evaluations=1)
    np.testing.assert_allclose(resumed.start, observer.events[-1].best_values)

    # Cancelled at a trial point whose calculation failed, the fit keeps the last one that worked.
    molecules = chargefw.MoleculeCollection((water(), formaldehyde()))
    template = make_parameter_template(molecules, "eem", atom_types="plain")
    plan = chargefw.assess(
        molecules, method="eem", parameter_set=template, execution="full"
    ).default_plan
    assert plan is not None
    failing = RecordingFitObserver(lambda progress: np.isnan(progress.objective))
    stopped = fit(
        FailingEvaluator(plan, threads=1),
        references,
        start=FIT_START,
        free=FIT_FREE,
        observer=failing,
    )
    assert np.isnan(failing.events[-1].objective) and failing.events[-1].failed_evaluations == 1
    assert not stopped.success
    np.testing.assert_array_equal(stopped.values, failing.events[-1].best_values)
    assert stopped.values[5] <= 2.9


def test_missing_references_leave_the_fit_and_the_errors() -> None:
    evaluator = training_evaluator()
    references = evaluator.charges(TRUE_VALUES)
    # The second water hydrogen has no reference.
    missing = references.copy()
    missing[2] = np.nan
    excluded = references.copy()
    excluded[2] = 5.0
    weights = np.ones(7)
    weights[2] = 0.0

    without = fit(evaluator, missing, start=FIT_START, free=FIT_FREE, max_evaluations=5)
    # A missing reference ignores its weight; a zero weight leaves a known reference out of the fit.
    weighted = fit(
        evaluator,
        missing,
        weights=np.where(np.isnan(missing), 7.0, 1.0),
        start=FIT_START,
        free=FIT_FREE,
        max_evaluations=5,
    )
    zero = fit(
        evaluator, excluded, weights=weights, start=FIT_START, free=FIT_FREE, max_evaluations=5
    )

    np.testing.assert_allclose(without.values, zero.values)
    np.testing.assert_allclose(without.values, weighted.values)
    charges = evaluator.charges(without.values)
    present = [0, 1, 3, 4, 5, 6]
    assert without.rmsd == pytest.approx(
        np.sqrt(np.mean((charges[present] - references[present]) ** 2))
    )
    assert without.diagnostics.by_element["H"].count == 3
    # A zero weight still reports the charge, against its reference.
    assert zero.diagnostics.by_element["H"].count == 4
    assert zero.rmsd > 1.0

    # Without any hydrogen reference, the hydrogen type leaves the errors but keeps its count; the
    # carbon charge, with weight 0 but a reference, is still reported.
    no_hydrogen = references.copy()
    no_hydrogen[[1, 2, 5, 6]] = np.nan
    without_carbon = np.ones(7)
    without_carbon[3] = 0.0
    partial = fit(
        evaluator,
        no_hydrogen,
        weights=without_carbon,
        start=FIT_START,
        free=FIT_FREE,
        max_evaluations=5,
    )
    assert {key.element for key in partial.diagnostics.by_atom_type} == {"C", "O"}
    assert set(partial.diagnostics.by_element) == {"C", "O"}
    assert partial.diagnostics.by_element["C"].count == 1
    assert partial.diagnostics.matches == without.diagnostics.matches

    infinite = references.copy()
    infinite[2] = np.inf
    with pytest.raises(ValueError, match="finite, or nan"):
        fit(evaluator, infinite, start=FIT_START, free=FIT_FREE)
    with pytest.raises(ValueError, match="no value at all"):
        fit(evaluator, np.full(7, np.nan), start=FIT_START, free=FIT_FREE)


def test_start_from_a_parameter_set_matches_bonds_in_either_order(tmp_path: Path) -> None:
    evaluator = template_evaluator("sqe", "hbo", "bo")
    references = evaluator.charges(SQE_TRUE_VALUES)
    written = fit(evaluator, references, start=SQE_TRUE_VALUES, free=[0], max_evaluations=1)
    path = tmp_path / "reversed.json"
    chargefw.write_parameter_set(written.parameter_set(id="reversed"), path)
    document = json.loads(path.read_text(encoding="utf-8"))
    for entry in document["bond"]["data"]:
        entry["key"]["atoms"].reverse()
    path.write_text(json.dumps(document), encoding="utf-8")

    template = template_evaluator("sqe", "hbo", "bo")
    started = fit(
        template,
        references,
        start=chargefw.read_parameter_set(path),
        free=[0],
        max_evaluations=1,
    )

    np.testing.assert_allclose(started.start, written.values)


def test_start_from_a_parameter_set_fills_only_the_slots_it_shares() -> None:
    evaluator = template_evaluator("sqe", "hbo", "bo")
    references = evaluator.charges(SQE_TRUE_VALUES)
    profile = PROFILES["sqe"].start

    def water_set(atom_types: str) -> chargefw.ParameterSet:
        # A set for water alone: none of formaldehyde's types.
        template = make_parameter_template(
            water(),
            "sqe",
            atom_types=atom_types,  # type: ignore[arg-type]
            bond_types="bo",
            id=f"water-{atom_types}",
        )
        plan = chargefw.assess(
            water(), method="sqe", parameter_set=template, execution="full"
        ).default_plan
        assert plan is not None
        water_evaluator = Evaluator(plan, threads=1)
        positions = {(slot.key, slot.name): index for index, slot in enumerate(evaluator.layout)}
        values = np.array(
            [
                SQE_TRUE_VALUES[positions[(slot.key, slot.name)]]
                if (slot.key, slot.name) in positions
                else profile[slot.name]
                for slot in water_evaluator.layout
            ]
        )
        written = fit(water_evaluator, water_evaluator.charges(values), start=values, free=[0])
        return written.parameter_set(id=f"water-{atom_types}-fitted")

    # The electronegativity of water's hydrogen is in both layouts, so the fit can vary it.
    free = [
        index
        for index, slot in enumerate(evaluator.layout)
        if slot.name == "electronegativity"
        and isinstance(slot.key, chargefw.AtomParameterKey)
        and slot.key.element == "H"
    ]
    hbo = water_set("hbo")
    shared = fit(evaluator, references, start=hbo, free=free, max_evaluations=1)
    water_keys = set(hbo.atom_keys) | set(hbo.bond_keys)
    for index, slot in enumerate(evaluator.layout):
        if slot.key in water_keys and index not in free:
            assert shared.start[index] == SQE_TRUE_VALUES[index], slot
        elif slot.key not in water_keys:
            assert shared.start[index] == profile[slot.name], slot

    # bonded types share no key with hbo types, so every slot starts as from the template.
    unrelated = fit(evaluator, references, start=water_set("bonded"), free=free, max_evaluations=1)
    np.testing.assert_array_equal(
        unrelated.start, [profile[slot.name] for slot in evaluator.layout]
    )


class FailingAfterCancellation(Evaluator):
    """Fails every calculation once ``failing`` is set."""

    failing = False

    def _normal_equations(
        self,
        values: NDArray[np.float64],
        references: NDArray[np.float64],
        weights: NDArray[np.float64],
    ) -> tuple[NDArray[np.float64], float, NDArray[np.float64], NDArray[np.float64]]:
        if self.failing:
            raise EvaluationError("calculation lost after cancellation")
        return super()._normal_equations(values, references, weights)


def test_cancelled_fit_whose_diagnostics_fail_leaves_the_checkpoint() -> None:
    molecules = chargefw.MoleculeCollection((water(), formaldehyde()))
    template = make_parameter_template(molecules, "eem", atom_types="plain")
    plan = chargefw.assess(
        molecules, method="eem", parameter_set=template, execution="full"
    ).default_plan
    assert plan is not None
    evaluator = FailingAfterCancellation(plan, threads=1)
    references = training_evaluator().charges(TRUE_VALUES)

    def stop_and_break(progress: FitProgress) -> bool:
        if progress.evaluation < 2:
            return False
        evaluator.failing = True
        return True

    observer = RecordingFitObserver(stop_and_break)
    with pytest.raises(EvaluationError, match="lost after cancellation"):
        fit(evaluator, references, start=FIT_START, free=FIT_FREE, observer=observer)

    checkpoint = observer.events[-1].parameter_set(id="checkpoint")
    resumed = fit(
        training_evaluator(), references, start=checkpoint, free=FIT_FREE, max_evaluations=1
    )
    np.testing.assert_allclose(resumed.start, observer.events[-1].best_values)


FIXED_ION_BETWEEN_WATERS = """data_ion_between_waters
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_formal_charge
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
HETATM 1 O O . HOH W . ? 0 0 0 1 20 0 5 HOH W O 1
HETATM 2 H H1 . HOH W . ? 0.96 0 0 1 20 0 5 HOH W H1 1
HETATM 3 H H2 . HOH W . ? -0.24 0.93 0 1 20 0 5 HOH W H2 1
HETATM 4 Mg MG . MG I . ? 0 0 3 1 20 0 9 MG I MG 1
HETATM 5 O O . HOH W . ? 0 0 6 1 20 0 6 HOH W O 1
HETATM 6 H H1 . HOH W . ? 0.96 0 6 1 20 0 6 HOH W H1 1
HETATM 7 H H2 . HOH W . ? -0.24 0.93 6 1 20 0 6 HOH W H2 1
#
"""


def test_atom_types_line_up_with_atoms_around_a_fixed_ion() -> None:
    molecules = chargefw.io.parse(FIXED_ION_BETWEEN_WATERS, format="mmcif", bonds="templates")
    assert [int(z) for z in molecules[0].atomic_numbers] == [8, 1, 1, 12, 8, 1, 1]
    template = make_parameter_template(molecules, "eem", atom_types="plain")
    plan = chargefw.assess(
        molecules, method="eem", parameter_set=template, execution="full", fixed_ions=["MG"]
    ).default_plan
    assert plan is not None
    evaluator = Evaluator(plan, threads=1)
    values = evaluator.values.copy()
    for index, slot in enumerate(evaluator.layout):
        values[index] = {"kappa": 0.5, "A": 2.0, "B": 1.0}[slot.name]
        if isinstance(slot.key, chargefw.AtomParameterKey) and slot.key.element == "O":
            values[index] = {"A": 3.0, "B": 1.0}[slot.name]
    charges = evaluator.charges(values)
    assert charges[3] == 2.0
    # Each element gets its own offset; misaligned types would mix them, or take in the ion.
    references = charges + np.array([0.01, 0.02, 0.02, 0.5, 0.01, 0.02, 0.02])

    result = fit(evaluator, references, start=values, free=[1], max_evaluations=1)

    by_type = {
        key.element: diagnostics for key, diagnostics in result.diagnostics.by_atom_type.items()
    }
    assert set(by_type) == {"H", "O"}
    assert (by_type["H"].count, by_type["O"].count) == (4, 2)
    fitted = evaluator.charges(result.values)
    hydrogens, oxygens = [1, 2, 5, 6], [0, 4]
    assert by_type["H"].mean_error == pytest.approx(
        np.mean(fitted[hydrogens] - references[hydrogens])
    )
    assert by_type["O"].mean_error == pytest.approx(np.mean(fitted[oxygens] - references[oxygens]))
