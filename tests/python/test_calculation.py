"""Focused calculation facade and result-model checks."""

import gc
import warnings
from collections.abc import Mapping
from concurrent.futures import ThreadPoolExecutor
from dataclasses import replace
from threading import Lock
from typing import Any, cast
from unittest.mock import patch

import chargefw
import numpy as np
import pytest
from chargefw.parameterization import Evaluator


def water(conformers: int = 1) -> chargefw.Molecule:
    coordinates = np.array(
        [
            [[0.0, 0.0, 0.0], [0.96, 0.0, 0.0], [-0.24, 0.93, 0.0]],
            [[0.0, 0.0, 0.1], [0.96, 0.0, 0.1], [-0.24, 0.93, 0.1]],
        ][:conformers]
    )
    return chargefw.Molecule(
        [8, 1, 1],
        bonds=np.array([[0, 1, 1], [0, 2, 1]], dtype=np.int32),
        coordinates=coordinates[0] if conformers == 1 else coordinates,
        source_name="water.sdf",
        record_index=2,
        record_id="water-record",
        atom_ids=["O", "H1", "H2"],
    )


FIXED_IONS_MMCIF = """data_fixed_groups
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
#
"""


def fixed_groups_molecules() -> chargefw.MoleculeCollection:
    return chargefw.io.parse(FIXED_IONS_MMCIF, format="mmcif", bonds="templates")


def calculate_formal(molecules: Any) -> chargefw.CalculationResult:
    return chargefw.calculate(molecules, method="formal", execution="full", threads=1)


def assess_formal(molecules: Any) -> chargefw.Assessment:
    return chargefw.assess(molecules, method="formal", execution="full", threads=1)


class RecordingObserver(chargefw.CalculationObserver):
    def __init__(self) -> None:
        self.events: list[chargefw.CalculationProgress] = []

    def on_progress(self, progress: chargefw.CalculationProgress) -> None:
        self.events.append(progress)


def test_all_ions_selects_noncommon_template_in_assessment_and_calculation() -> None:
    # CU1's atom name is CU and its template charge differs from the imported zero.
    text = FIXED_IONS_MMCIF.replace(
        "HETATM 4 Mg MG . MG I . ? 0 0 3 1 20 0 9 MG I MG 1",
        "HETATM 4 Cu CU . CU1 I . ? 0 0 3 1 20 0 9 CU1 I CU 1",
    )
    molecules = chargefw.io.parse(text, format="mmcif", bonds="templates")
    common = chargefw.calculate(molecules, method="formal", fixed_ions=chargefw.COMMON_IONS)
    assert common.plan is not None, (
        "unselected noncommon ion must leave formal calculation available"
    )
    assert common.plan.fixed_ions is None
    np.testing.assert_array_equal(common.assignments[0].values, [0.0] * 4)
    assessment = chargefw.assess(
        molecules,
        method="sqeqp",
        parameter_set="SQEqp_Schindler2021_CCD_gen",
        execution="full",
        fixed_ions=chargefw.ALL_IONS,
    )
    assert assessment.default_plan is not None, (
        "ALL_IONS must produce an SQEqp plan with the CU1 template"
    )
    result = assessment.default_plan.calculate()
    assert result.requested.fixed_ions == chargefw.ALL_IONS
    assert result.assignments[0].values[3] == 1.0
    assert float(result.assignments[0].values[:3].sum()) == pytest.approx(0.0, abs=1e-12)
    np.testing.assert_array_equal(molecules[0].formal_charges, [0] * 4)
    assert result.plan is not None and result.plan.fixed_ions is not None, (
        "ALL_IONS must retain effective source provenance"
    )
    assert result.plan.fixed_ions.sources == (chargefw.FixedAtomCharge(0, 3, 1.0, "CU1"),)


def test_observer_receives_owned_execution_progress() -> None:
    observer = RecordingObserver()
    collection = chargefw.MoleculeCollection([chargefw.Molecule([1]), chargefw.Molecule([8, 1])])
    result = chargefw.calculate(
        collection,
        method="formal",
        execution="full",
        threads=1,
        observer=observer,
    )

    assert result.status == "success"
    assert [event.phase for event in observer.events] == [
        "computation_started",
        "target_started",
        "target_finished",
        "target_started",
        "target_finished",
        "computation_finished",
    ]
    started = observer.events[0]
    assert started.method.id == "formal"
    assert started.mode == "full"
    targets = [event for event in observer.events if event.phase == "target_started"]
    assert [event.target_index for event in targets] == [0, 1]
    assert [event.target_count for event in targets] == [2, 2]
    assert [event.molecule_index for event in targets] == [0, 1]
    assert [event.conformer_index for event in targets] == [None, None]
    assert all(event.elapsed_seconds >= 0.0 for event in observer.events)


def test_reduced_execution_reports_policy_and_fragment_progress() -> None:
    for mode in ("cutoff", "cover"):
        observer = RecordingObserver()
        result = chargefw.calculate(
            water(),
            method="eem",
            execution=cast(Any, mode),
            radius=8.0,
            threads=1,
            observer=observer,
        )

        assert result.status == "success"
        assert result.plan is not None, "successful calculation must report plan provenance"
        assert result.plan.execution == chargefw.Execution(mode=mode, radius=8.0)
        fragments = [event for event in observer.events if event.phase == "fragment_progress"]
        assert fragments
        assert fragments[-1].completed_fragment_count == fragments[-1].fragment_count
        assert fragments[-1].target_index == 0
        assert fragments[-1].target_count == 1


def test_fixed_ions_snapshot_resolve_and_reuse_plan() -> None:
    molecules = fixed_groups_molecules()
    assert molecules[0].bond_count == 2
    names = ["MG"]
    assessment = chargefw.assess(
        molecules,
        method="sqeqp",
        parameter_set="SQEqp_Schindler2021_CCD_gen",
        execution="full",
        fixed_ions=names,
    )
    names.append("CA")
    plan = assessment.default_plan
    assert plan is not None, "assessment must produce an SQEqp plan with the selected MG source"
    del assessment
    gc.collect()
    result = plan.calculate()
    direct = chargefw.calculate(
        molecules,
        method="sqeqp",
        parameter_set="SQEqp_Schindler2021_CCD_gen",
        execution="full",
        fixed_ions=chargefw.COMMON_IONS,
    )
    assert result.requested.fixed_ions == ("MG",)
    assert result.assignments[0].values[3] == 2.0
    np.testing.assert_allclose(result.assignments[0].values, direct.assignments[0].values)
    assert molecules[0].formal_charges.tolist() == [0, 0, 0, 0]
    assert result.plan is not None and result.plan.fixed_ions is not None, (
        "executed plan must retain effective source provenance"
    )
    assert result.plan.fixed_ions.sources == (chargefw.FixedAtomCharge(0, 3, 2.0, "MG"),)
    assert result.plan.fixed_ions == chargefw.FixedIons([chargefw.FixedAtomCharge(0, 3, 2.0, "MG")])
    with pytest.raises(AttributeError):
        cast(Any, result.plan.fixed_ions).sources = ()
    with pytest.raises(AttributeError):
        cast(Any, result.requested).fixed_ions = ("CA",)


def test_fixed_ions_validation_empty_and_method_support() -> None:
    molecules = fixed_groups_molecules()
    for invalid in ("MG", b"MG", ["MG", 12]):
        with pytest.raises(TypeError):
            chargefw.assess(molecules, fixed_ions=cast(Any, invalid))
    with pytest.raises(ValueError):
        chargefw.assess(molecules, fixed_ions=["mg"])
    with pytest.raises(ValueError):
        chargefw.assess(chargefw.MoleculeCollection([]), fixed_ions=["bad"])

    absent = chargefw.calculate(molecules, method="formal", fixed_ions=["CA"])
    assert absent.requested.fixed_ions == ("CA",)
    assert absent.plan is not None, "formal calculation must produce a plan"
    assert absent.plan.fixed_ions is None
    duplicate = chargefw.calculate(molecules, method="eem", fixed_ions=["MG", "MG"])
    assert duplicate.plan is not None and duplicate.plan.fixed_ions is not None, (
        "duplicate names must resolve to one fixed source"
    )
    assert len(duplicate.plan.fixed_ions.sources) == 1
    assert chargefw.RequestedCalculation().fixed_ions == ()

    plain = chargefw.calculate(water(), method="formal")
    empty = chargefw.calculate(water(), method="formal", fixed_ions=[])
    assert empty.requested.fixed_ions == ()
    np.testing.assert_array_equal(plain.assignments[0].values, empty.assignments[0].values)

    unsupported = chargefw.assess(
        molecules, method="qeq", parameter_set="QEq_original", fixed_ions=["MG"]
    )
    assert not unsupported.plans
    assert "unsupported_fixed_ions" in [
        issue.kind for rejection in unsupported.rejections for issue in rejection.issues
    ]
    with pytest.raises(ValueError):
        chargefw.assess(water(), fixed_ions=["MG"])


def test_observer_receives_parallel_target_progress() -> None:
    class CountingObserver(chargefw.CalculationObserver):
        def __init__(self) -> None:
            self.lock = Lock()
            self.completed = 0

        def on_progress(self, progress: chargefw.CalculationProgress) -> None:
            if progress.phase == "target_finished":
                with self.lock:
                    self.completed += 1

    molecule_count = 64
    observer = CountingObserver()
    result = chargefw.calculate(
        chargefw.MoleculeCollection(chargefw.Molecule([1]) for _ in range(molecule_count)),
        method="formal",
        execution="full",
        threads=2,
        observer=observer,
    )

    assert result.status == "success"
    assert observer.completed == molecule_count


def test_observer_can_cancel_a_reusable_plan() -> None:
    class CancellingObserver(chargefw.CalculationObserver):
        def __init__(self) -> None:
            self.cancel_requested = False
            self.phases: list[str] = []

        def on_progress(self, progress: chargefw.CalculationProgress) -> None:
            self.phases.append(progress.phase)
            if progress.phase == "target_started":
                self.cancel_requested = True

        def cancelled(self) -> bool:
            return self.cancel_requested

    cases = (
        (chargefw.MoleculeCollection([water()]), None),
        (fixed_groups_molecules(), ["MG"]),
    )
    for molecule, groups in cases:
        observer = CancellingObserver()
        options: dict[str, Any] = {
            "method": "eem" if groups is None else "sqeqp",
            "execution": "full",
            "fixed_ions": groups,
        }
        if groups is not None:
            options["parameter_set"] = "SQEqp_Schindler2021_CCD_gen"
        assessment = chargefw.assess(molecule, **options)
        plan = assessment.default_plan
        assert plan is not None, "assessment must produce a default plan"
        with pytest.raises(chargefw.CalculationCancelledError) as raised:
            plan.calculate(threads=1, observer=observer)

        cancelled = raised.value.result
        assert cancelled.status == "cancelled"
        assert cancelled.assignments == ()
        assert observer.phases[0] == "computation_started"
        assert "target_started" in observer.phases
        assert observer.phases[-1] == "computation_finished"
        assert cancelled.requested.fixed_ions == tuple(groups or ())
        if groups is not None:
            assert cancelled.plan is not None
            assert cancelled.plan.fixed_ions is not None
            assert cancelled.plan.fixed_ions.sources == (chargefw.FixedAtomCharge(0, 3, 2.0, "MG"),)

        result = plan.calculate(threads=1)
        assert result.status == "success"
        assert len(result.assignments) == molecule[0].conformer_count
        assert all(
            assignment.values.shape == (molecule[0].atom_count,)
            for assignment in result.assignments
        )
        assert all(np.all(np.isfinite(assignment.values)) for assignment in result.assignments)


def test_observer_must_use_public_base_class() -> None:
    with pytest.raises(TypeError, match="CalculationObserver"):
        chargefw.calculate(water(), observer=cast(Any, object()))


def test_observer_callback_failure_does_not_change_calculation() -> None:
    class FailingObserver(chargefw.CalculationObserver):
        def on_progress(self, progress: chargefw.CalculationProgress) -> None:
            raise RuntimeError(f"failed during {progress.phase}")

    with patch("sys.unraisablehook") as unraisable:
        result = chargefw.calculate(
            water(),
            method="formal",
            execution="full",
            threads=1,
            observer=FailingObserver(),
        )

    assert result.status == "success"
    assert unraisable.called
    # The recorded exception owns its traceback, which otherwise retains this calculation frame
    # until interpreter shutdown and produces a false-positive Nanobind leak report.
    unraisable.reset_mock()


def test_method_scope_controls_methods_that_are_not_recommended() -> None:
    assert not chargefw.methods["veem"].recommended
    assert chargefw.methods["eem"].recommended

    # VEEM needs no parameters and applies to water, but it is not recommended.
    recommended = chargefw.assess(water())
    assert "veem" not in {plan.method.id for plan in recommended.plans}

    every = chargefw.assess(water(), method_scope="all")
    assert "veem" in {plan.method.id for plan in every.plans}
    assert chargefw.calculate(water(), method_scope="all").requested.method_scope == "all"
    default = chargefw.calculate(water())
    assert default.status == "success"
    assert default.plan is not None
    assert default.requested.method_scope == "recommended"

    named = chargefw.calculate(water(), method="veem")
    assert named.plan is not None and named.plan.method.id == "veem"

    with pytest.raises(ValueError):
        chargefw.assess(water(), method="veem", method_scope="all")
    with pytest.raises(ValueError):
        chargefw.assess(water(), method_scope=cast(Any, "some"))


def test_assessment_and_full_calculation() -> None:
    molecule = water()
    assessment = chargefw.assess(molecule, method="eem", execution="full")
    assert assessment.plans
    assert all(plan.execution.mode == "full" for plan in assessment.plans)
    plan = assessment.default_plan
    assert plan is not None, "assessment must produce a default plan"
    assert plan.execution == chargefw.Execution(
        mode="full",
        radius=None,
    )
    assert plan.method.id == "eem"

    result = plan.calculate()
    assert result.status == "success"
    assert len(result.assignments) == 1
    assignment = result.assignments[0]
    assert assignment.molecule_index == 0
    assert assignment.conformer_index == 0
    assert assignment.source.record_id == "water-record"
    assert assignment.atom_ids == ("O", "H1", "H2")
    assert assignment.values.dtype == np.dtype(np.float64)
    assert assignment.values.flags.c_contiguous
    assert not assignment.values.flags.writeable
    assert assignment.values.ctypes.data % assignment.values.dtype.alignment == 0
    with pytest.raises(ValueError):
        assignment.values.setflags(write=True)
    assert np.isclose(assignment.values.sum(), 0.0)
    assert np.asarray(assignment) is assignment.values
    assert np.asarray(assignment, dtype=np.float32).dtype == np.dtype(np.float32)
    assert result.requested.method is chargefw.methods["eem"]
    assert result.requested.execution == "full"
    executed = result.plan
    assert executed is not None, "successful calculation must report plan provenance"
    assert executed.method.id == "eem"
    assert executed.parameter_set == plan.parameter_set
    assert executed.execution.mode == "full"
    assert result.timings.applicability_seconds >= 0.0
    assert result.timings.computation_seconds >= 0.0
    with pytest.raises(TypeError):
        cast(Any, result.requested.method_options)["iters"] = 1

    overridden_threads = plan.calculate(threads=1)
    assert overridden_threads.requested.threads == 1


def test_plans_are_filtered_and_reusable() -> None:
    molecule = water()
    assessment = chargefw.assess(molecule, method="eqeq")
    assert {candidate_plan.execution.mode for candidate_plan in assessment.plans} == {
        "full",
        "cutoff",
        "cover",
    }
    for candidate_plan in assessment.plans:
        assert candidate_plan.calculate(threads=1).status == "success"

    cutoff = chargefw.assess(
        molecule,
        method="eqeq",
        execution="cutoff",
        radius=8.0,
    )
    assert cutoff.plans
    assert all(plan.execution.mode == "cutoff" for plan in cutoff.plans)

    default_plan = assessment.default_plan
    assert default_plan is not None, "assessment must produce a default plan"
    with ThreadPoolExecutor(max_workers=3) as executor:
        results = list(
            executor.map(
                lambda _: default_plan.calculate(threads=1),
                range(3),
            )
        )
    assert all(result.status == "success" for result in results)


def test_settings_apply_exactly_the_choices_of_a_plan_to_other_molecules() -> None:
    # The threshold makes automatic execution choose cutoff; CA is absent from the molecules.
    reference = chargefw.calculate(
        fixed_groups_molecules(),
        method="eem",
        cutoff_threshold=1,
        parameter_matching="permissive",
        fixed_ions=["MG", "CA"],
    )
    assert reference.plan is not None
    settings = reference.plan.settings
    assert settings.execution == chargefw.Execution(mode="cutoff", radius=12.0)
    assert settings.fixed_ions == ("MG", "CA")

    # With the default thresholds, automatic execution would choose full.
    plan = chargefw.assess(fixed_groups_molecules(), settings=settings).default_plan
    assert plan is not None
    assert plan.settings == settings
    result = plan.calculate(threads=1)
    assert result.plan is not None and result.plan.fixed_ions is not None
    assert result.plan.settings == settings
    assert [source.component_id for source in result.plan.fixed_ions.sources] == ["MG"]

    # Also an argument that has its default value is a choice that settings cannot be combined with.
    conflicts: list[dict[str, Any]] = [
        {"execution": "full"},
        {"execution": "auto"},
        {"parameter_matching": "strict"},
        {"method_scope": "recommended"},
        {"fixed_ions": []},
    ]
    for function in (chargefw.calculate, chargefw.assess):
        for conflict in conflicts:
            with pytest.raises(ValueError):
                cast(Any, function)(fixed_groups_molecules(), settings=settings, **conflict)

    # Automatic selection finds a method for xenon; the settings of an EEM plan do not fall back.
    xenon = chargefw.Molecule([54], coordinates=[[0.0, 0.0, 0.0]])
    water_plan = chargefw.calculate(water(), method="eem").plan
    assert water_plan is not None
    assert chargefw.calculate(xenon).status == "success"
    with pytest.raises(chargefw.NoExecutablePlanError) as context:
        chargefw.calculate(xenon, settings=water_plan.settings)
    assert [
        (rejection.method, rejection.parameter_set) for rejection in context.value.result.rejections
    ] == [(water_plan.method, water_plan.parameter_set)]


def test_objects_that_chargefw_creates_cannot_be_constructed() -> None:
    # Arguments that the constructors accepted when they were public.
    bundled = chargefw.parameter_sets["EEM_Cheminf_b3lyp_npa"]
    with pytest.raises(TypeError):
        cast(Any, chargefw.ParameterSet)(
            id=bundled.id,
            method=bundled.method.id,
            name=bundled.name,
            publication=bundled.publication,
            notes=bundled.notes,
            atom_keys=bundled.atom_keys,
            bond_keys=bundled.bond_keys,
        )
    result = chargefw.calculate(water(), method="eem")
    assert result.plan is not None
    settings = result.plan.settings
    with pytest.raises(TypeError):
        cast(Any, chargefw.CalculationSettings)(
            settings.method,
            settings.parameter_set,
            settings.method_options,
            settings.execution,
            result.requested,
        )
    # A copy with another ID would name a set that does not exist.
    with pytest.raises(TypeError):
        replace(cast(Any, bundled), id="mine")


def test_plan_estimates_resources_for_a_thread_limit() -> None:
    molecules = [water(), water()]
    plan = chargefw.assess(molecules, method="eqeq", execution="full").default_plan
    assert plan is not None, "assessment must produce a default plan"

    serial = plan.estimate(threads=1)
    parallel = plan.estimate(threads=2)
    assert serial.memory_bytes > 0.0
    assert serial.operations > 0.0
    assert parallel.memory_bytes > serial.memory_bytes
    assert parallel.operations == serial.operations
    assert plan.estimate().memory_bytes >= serial.memory_bytes
    with pytest.raises(ValueError):
        plan.estimate(threads=-1)


def test_assignment_cardinality_and_mapping() -> None:
    geometry_independent = chargefw.calculate(
        chargefw.Molecule([8, 1, 1]),
        method="formal",
        execution="full",
    )
    assert geometry_independent.status == "success"
    assert len(geometry_independent.assignments) == 1
    assert geometry_independent.assignments[0].conformer_index is None
    assert geometry_independent.assignment(molecule=0) is geometry_independent.assignments[0]
    with pytest.raises(ValueError):
        geometry_independent.assignment(molecule=0, conformer=0)

    collection = chargefw.MoleculeCollection([water(2)])
    multi_result = chargefw.calculate(collection, method="qeq", execution="full")
    assert multi_result.status == "success"
    assert [item.conformer_index for item in multi_result.assignments] == [0, 1]
    assert multi_result.assignments_by_molecule == (multi_result.assignments,)
    assert multi_result.assignment(molecule=0, conformer=0) is multi_result.assignments[0]
    assert multi_result.assignment(molecule=0, conformer=1) is multi_result.assignments[1]
    invalid_indices: tuple[Any, ...] = (True, np.bool_(False))
    for value in invalid_indices:
        with pytest.raises(TypeError):
            multi_result.assignment(molecule=value, conformer=0)
        with pytest.raises(TypeError):
            multi_result.assignment(molecule=0, conformer=value)
    assert (
        multi_result.assignment(molecule=cast(Any, np.int64(0)), conformer=cast(Any, np.int64(1)))
        is multi_result.assignments[1]
    )
    with pytest.raises(ValueError):
        multi_result.assignment(molecule=0)

    mapped_collection = chargefw.MoleculeCollection(
        [
            chargefw.Molecule([1], source_name="first", record_id="record-a", atom_ids=["A"]),
            chargefw.Molecule(
                [8, 1],
                source_name="second",
                record_id="record-b",
                atom_ids=["B", "C"],
            ),
        ]
    )
    mapped_result = calculate_formal(mapped_collection)
    assert mapped_result.status == "success"
    assert [assignment.molecule_index for assignment in mapped_result.assignments] == [0, 1]
    assert [assignment.source.record_id for assignment in mapped_result.assignments] == [
        "record-a",
        "record-b",
    ]
    assert [assignment.atom_ids for assignment in mapped_result.assignments] == [("A",), ("B", "C")]
    assert mapped_result.assignments_by_molecule == (
        (mapped_result.assignments[0],),
        (mapped_result.assignments[1],),
    )


def test_public_charge_assignment_is_validated_and_comparable() -> None:
    source = chargefw.SourceIdentity(source="fixture", record_index=1, record_id="record")
    assignment = chargefw.ChargeAssignment(
        values=np.array([0.25, -0.25]),
        molecule_index=cast(Any, np.int64(2)),
        conformer_index=0,
        source=source,
        atom_ids=("A", "B"),
    )
    same = chargefw.ChargeAssignment(
        values=np.array([0.25, -0.25]),
        molecule_index=2,
        conformer_index=0,
        source=source,
        atom_ids=("A", "B"),
    )
    different = chargefw.ChargeAssignment(
        values=np.array([0.5, -0.5]),
        molecule_index=2,
        conformer_index=0,
        source=source,
        atom_ids=("A", "B"),
    )
    assert assignment == same
    assert assignment != different
    assert not assignment.values.flags.writeable
    with pytest.raises(ValueError):
        assignment.values.setflags(write=True)

    invalid_assignments = (
        (ValueError, {"values": [[0.0]], "atom_ids": ("A",)}),
        (ValueError, {"values": [np.nan], "atom_ids": ("A",)}),
        (ValueError, {"values": [0.0], "atom_ids": ()}),
        (
            ValueError,
            {"values": [0.0], "atom_ids": ("A",), "molecule_index": -1},
        ),
        (
            TypeError,
            {"values": [0.0], "atom_ids": ("A",), "source": "fixture"},
        ),
    )
    defaults: dict[str, Any] = {
        "molecule_index": 0,
        "conformer_index": None,
        "source": source,
    }
    for error_type, overrides in invalid_assignments:
        with pytest.raises(error_type):
            chargefw.ChargeAssignment(**(defaults | overrides))


def test_functional_api_supports_concurrent_calculations() -> None:
    collection = chargefw.MoleculeCollection([chargefw.Molecule([8, 1, 1])])
    with ThreadPoolExecutor(max_workers=4) as executor:
        futures = [executor.submit(calculate_formal, collection) for _ in range(4)]
    results = [future.result() for future in futures]
    assert all(result.status == "success" for result in results)
    assert [result.assignments[0].values.tolist() for result in results] == [[0.0, 0.0, 0.0]] * 4


def test_automatic_thresholds_and_explicit_full_warnings() -> None:
    automatic = chargefw.assess(
        water(),
        method="eem",
        cutoff_threshold=1,
        cover_threshold=100,
    )
    assert automatic.default_plan is not None, "assessment must produce a default plan"
    assert automatic.default_plan.execution.mode == "cutoff"
    assert any(
        rejection.execution is not None and rejection.execution.mode == "full"
        for rejection in automatic.rejections
    )

    explicit_full = chargefw.assess(
        water(),
        method="eem",
        execution="full",
        cutoff_threshold=1,
        cover_threshold=100,
    )
    assert explicit_full.default_plan is not None, "explicit full assessment must produce a plan"
    assert [warning.kind for warning in explicit_full.default_plan.warnings] == [
        "resource_threshold_exceeded"
    ]


def test_molecules_without_hydrogens_are_reported() -> None:
    # A molecule built without hydrogens, as from RDKit without AddHs(), next to water.
    molecules = [chargefw.Molecule([6, 8], bonds=[[0, 1, 3]]), water()]
    assessment = chargefw.assess(molecules, method="formal")
    with pytest.warns(chargefw.InputWarning, match=r"^molecule 1: .*\(hydrogens_missing; "):
        result = chargefw.calculate(molecules, method="formal")

    for issues in (assessment.input_warnings, result.input_warnings):
        assert [(issue.kind, issue.molecule_index, issue.conformer_index) for issue in issues] == [
            ("hydrogens_missing", 0, None)
        ]


def test_close_atoms_are_reported_with_their_atom() -> None:
    # The first hydrogen is 0.3 Å from the oxygen; EQeq uses coordinates.
    molecule = chargefw.Molecule(
        [8, 1, 1],
        bonds=[[0, 1, 1], [0, 2, 1]],
        coordinates=[[0.0, 0.0, 0.0], [0.3, 0.0, 0.0], [-0.24, 0.93, 0.0]],
    )
    assert chargefw.assess(molecule, method="eqeq").input_warnings == ()
    with pytest.warns(chargefw.InputWarning, match=r"\(atoms_too_close; "):
        result = chargefw.calculate(molecule, method="eqeq")

    assert [
        (issue.kind, issue.molecule_index, issue.conformer_index, issue.atom_index)
        for issue in result.input_warnings
    ] == [("atoms_too_close", 0, 0, 1)]
    assert result.input_warnings[0].atom_indices.tolist() == [0, 1]


def test_issues_that_can_make_charges_wrong_raise_python_warnings() -> None:
    def chargefw_warnings(call: Any) -> list[warnings.WarningMessage]:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            call()
        return [item for item in caught if issubclass(item.category, chargefw.ChargeFWWarning)]

    # Only issues that can make the charges wrong warn; reading warns where the caller reads.
    (unbonded,) = chargefw_warnings(
        lambda: chargefw.io.parse(FIXED_IONS_MMCIF, format="mmcif", bonds="explicit")
    )
    assert unbonded.category is chargefw.InputWarning
    assert "(unbonded_residue_atoms; see molecule.import_issues)" in str(unbonded.message)
    assert unbonded.filename == __file__
    assert chargefw_warnings(lambda: calculate_formal(water())) == []
    assert chargefw_warnings(lambda: chargefw.assess(chargefw.Molecule([6]))) == []

    # A batch gets one warning per kind, which counts the molecules.
    carbon_monoxide = chargefw.Molecule([6, 8], bonds=[[0, 1, 3]])
    (batch,) = chargefw_warnings(
        lambda: calculate_formal([carbon_monoxide, water(), carbon_monoxide])
    )
    assert str(batch.message).startswith("2 of 3 molecules: hydrogens_missing, for example")

    # An assessed plan warns where the caller calculates it.
    plan = assess_formal(carbon_monoxide).default_plan
    assert plan is not None
    (from_plan,) = chargefw_warnings(plan.calculate)
    assert from_plan.category is chargefw.InputWarning
    assert from_plan.filename == __file__

    # A calculation that fails raises its own error, also when warnings are errors; the issues are
    # in its result.
    with warnings.catch_warnings():
        warnings.simplefilter("error", chargefw.ChargeFWWarning)
        with pytest.raises(chargefw.NoExecutablePlanError) as context:
            chargefw.calculate(carbon_monoxide, method="qeq", execution="full")
    assert [issue.kind for issue in context.value.result.input_warnings] == ["hydrogens_missing"]


def test_no_plan_result_and_typed_exception() -> None:
    molecule = chargefw.Molecule([8, 1, 1], bonds=[[0, 1, 1], [0, 2, 1]])
    assessment = chargefw.assess(
        molecule,
        method="qeq",
        execution="full",
    )
    assert not assessment.plans
    with pytest.raises(chargefw.NoExecutablePlanError) as context:
        chargefw.calculate(molecule, method="qeq", execution="full")
    result = context.value.result
    assert result.status == "no_executable_plan"
    assert result.assignments == ()
    assert result.rejections
    assert result.rejections[0].issues[0].kind == "missing_feature"
    assert result.rejections[0].issues[0].message in str(context.value)


def test_nonfinite_coordinates_are_assessed_per_method() -> None:
    molecule = chargefw.Molecule(
        [8, 1, 1],
        bonds=[[0, 1, 1], [0, 2, 1]],
        coordinates=[[0.0, 0.0, 0.0], [np.nan, 0.0, 0.0], [-0.24, 0.93, 0.0]],
    )
    assert np.isnan(molecule.coordinates[0, 1, 0])

    formal_result = chargefw.calculate(molecule, method="formal", execution="full")
    assert formal_result.status == "success"

    assessment = chargefw.assess(
        molecule,
        method="qeq",
        parameter_set="QEq_original",
        execution="full",
    )
    assert not assessment.plans
    invalid_geometry = [
        issue
        for rejection in assessment.rejections
        for issue in rejection.issues
        if issue.kind == "invalid_geometry"
    ]
    assert len(invalid_geometry) == 1
    assert invalid_geometry[0].atom_index == 1
    assert invalid_geometry[0].conformer_index == 0


def test_invalid_selection_requests_raise_value_error() -> None:
    with pytest.raises(ValueError):
        chargefw.assess(water(), method="not-a-method")

    with pytest.raises(ValueError):
        chargefw.assess(water(), parameter_set="QEq_original")

    with pytest.raises(ValueError):
        chargefw.assess(water(), method="eem", parameter_set="not-a-parameter-set")

    with pytest.raises(ValueError):
        chargefw.assess(water(), method="eem", parameter_set="QEq_original")


def test_results_outlive_calculation_inputs() -> None:
    def calculate_owned_result() -> tuple[chargefw.CalculationResult, chargefw.Plan]:
        molecule = chargefw.Molecule(
            [8, 1, 1],
            source_name="owned",
            record_id="owned-record",
            atom_ids=["O", "H1", "H2"],
        )
        assessment = assess_formal(molecule)
        plan = assessment.default_plan
        assert plan is not None, "assessment must produce a default plan"
        return plan.calculate(), plan

    result, plan = calculate_owned_result()
    gc.collect()
    assert plan.execution.mode == "full"
    assert result.status == "success"
    assert result.assignments[0].source.record_id == "owned-record"
    assert result.assignments[0].atom_ids == ("O", "H1", "H2")
    np.testing.assert_array_equal(result.assignments[0].values, [0.0, 0.0, 0.0])

    values = result.assignments[0].values
    del result
    gc.collect()
    np.testing.assert_array_equal(values, [0.0, 0.0, 0.0])
    assert not values.flags.writeable


def test_invalid_options_are_rejected_early() -> None:
    invalid_options = (
        (ValueError, {"execution": "fast"}),
        (ValueError, {"parameter_matching": "guess"}),
        (TypeError, {"threads": True}),
        (TypeError, {"cutoff_threshold": False}),
        (ValueError, {"cover_threshold": np.iinfo(np.uintp).max + 1}),
    )
    for error_type, options in invalid_options:
        with pytest.raises(error_type):
            chargefw.assess(water(), **options)

    # Options belong to a named method.
    with pytest.raises(ValueError):
        chargefw.assess(water(), method_options={"iters": 2})
    # Also an empty value that is not a mapping is refused, not taken as no options.
    not_mappings: tuple[object, ...] = ([], "", 0, False)
    for not_a_mapping in not_mappings:
        with pytest.raises(TypeError):
            chargefw.assess(water(), method="peoe", method_options=cast(Any, not_a_mapping))


def test_method_option_overrides_are_validated_and_reported() -> None:
    molecule = chargefw.Molecule([8, 1, 1], bonds=[[0, 1, 1], [0, 2, 1]])
    result = chargefw.calculate(
        molecule,
        method="peoe",
        method_options={"iters": 2},
        execution="full",
    )
    assert result.status == "success"
    assert result.plan is not None, "successful calculation must report plan provenance"
    assert dict(result.requested.method_options) == {"iters": 2}
    assert dict(result.plan.method_options) == {"initial_charges": "zero", "iters": 2}

    for options in (
        {"iters": 0},
        {"iters": 26},
        {"initial_charges": "unknown"},
        {"unknown": 1},
    ):
        with pytest.raises(ValueError):
            chargefw.assess(
                molecule,
                method="peoe",
                method_options=options,
                execution="full",
            )


def test_catalogs_and_descriptors_are_immutable_values() -> None:
    methods = chargefw.methods
    assert isinstance(methods, Mapping)
    assert "eem" in methods
    assert methods.get("eem") == methods["eem"]
    assert methods.get("not-a-method") is None
    invalid_keys: tuple[object, ...] = (5, None, [])
    for catalog in (methods, chargefw.parameter_sets):
        assert "not-an-id" not in catalog
        for key in invalid_keys:
            assert key not in catalog
            with pytest.raises(TypeError):
                catalog[key]  # type: ignore[index]
    assert tuple(methods) == tuple(method.id for method in methods.values())
    assert methods["eem"].id == "eem"
    with pytest.raises(KeyError):
        methods["not-a-method"]
    with pytest.raises(TypeError):
        cast(Any, methods)["unexpected"] = methods["eem"]
    with pytest.raises(TypeError):
        del cast(Any, methods)["eem"]

    eem = methods["eem"]
    peoe = methods["peoe"]
    assert peoe.requires_bonds and not eem.requires_bonds
    assert isinstance(peoe.options, Mapping)
    assert peoe.options["iters"].id == "iters"

    assert isinstance(chargefw.parameter_sets, Mapping)
    assert eem.parameter_sets
    assert all(parameter_set.method.id == "eem" for parameter_set in eem.parameter_sets.values())
    parameter_set = next(iter(eem.parameter_sets.values()))
    assert chargefw.parameter_sets[parameter_set.id] == parameter_set

    qeq = chargefw.parameter_sets["QEq_original"]
    assert {key.classifier for key in qeq.atom_keys} == {"plain"}
    assert qeq.bond_keys == ()
    peptide = chargefw.parameter_sets["SQEqp_Schindler2021_PUB_pept"]
    assert (
        chargefw.AtomParameterKey(element="C", classifier="bonded", type="COO") in peptide.atom_keys
    )
    bond = peptide.bond_keys[0]
    assert bond.classifier == "bo"
    assert all(atom.classifier == "bonded" for atom in bond.atoms)
    result = chargefw.calculate(
        water(),
        method=eem,
        parameter_set=parameter_set,
        execution="full",
    )
    assert result.status == "success"
    executed = result.plan
    assert executed is not None, "successful calculation must report plan provenance"
    assert executed.parameter_set == parameter_set


def write_eem_set(path: Any, set_id: str, oxygen_a: float) -> Any:
    path.write_text(
        f"""{{
  "metadata": {{"id": "{set_id}", "method": "eem", "name": "Water EEM"}},
  "common": {{"names": ["kappa"], "values": [0.5]}},
  "atom": {{
    "names": ["A", "B"],
    "data": [
      {{"key": {{"element": "H", "classifier": "plain", "type": "*"}}, "values": [2.0, 1.0]}},
      {{"key": {{"element": "O", "classifier": "plain", "type": "*"}}, "values": [{oxygen_a}, 1.0]}}
    ]
  }}
}}
"""
    )
    return path


def test_parameter_set_read_from_file_is_used_and_reported(tmp_path: Any) -> None:
    equal = chargefw.read_parameter_set(write_eem_set(tmp_path / "equal.json", "water-equal", 2.0))
    polar = chargefw.read_parameter_set(write_eem_set(tmp_path / "polar.json", "water-polar", 3.0))
    assert equal.id == "water-equal" and equal.method.id == "eem"
    assert "water-equal" not in chargefw.parameter_sets

    # Equal electronegativities give zero charges; a higher one for oxygen makes it negative.
    neutral = chargefw.calculate(water(), method="eem", parameter_set=equal, threads=1)
    np.testing.assert_allclose(neutral.assignments[0].values, [0.0] * 3, atol=1e-12)
    charged = chargefw.calculate(water(), method="eem", parameter_set=polar, threads=1)
    assert charged.assignments[0].values[0] < 0.0
    assert charged.plan is not None and charged.plan.parameter_set == polar

    molecule = water()
    assessment = chargefw.assess(molecule, method="eem", parameter_set=polar, threads=1)
    plan = assessment.default_plan
    assert plan is not None and plan.parameter_set == polar
    from_plan = plan.calculate(threads=1)
    assert from_plan.plan is not None and from_plan.plan.parameter_set == polar
    np.testing.assert_array_equal(from_plan.assignments[0].values, charged.assignments[0].values)


def test_parameter_sets_read_from_files_with_other_values_are_distinct(tmp_path: Any) -> None:
    first = chargefw.read_parameter_set(write_eem_set(tmp_path / "first.json", "water", 2.0))
    second = chargefw.read_parameter_set(write_eem_set(tmp_path / "second.json", "water", 3.0))

    # Same ID, metadata, and keys; only the values differ.
    assert (first.id, first.atom_keys) == (second.id, second.atom_keys)
    assert first != second
    assert first == first
    assert len({first, second}) == 2
    assert chargefw.parameter_sets["EEM_Baek1991"] == chargefw.parameter_sets["EEM_Baek1991"]


def test_parameter_set_read_from_file_rejects_bundled_ids_and_bad_files(tmp_path: Any) -> None:
    bundled_id = next(iter(chargefw.parameter_sets))
    with pytest.raises(ValueError, match="bundled"):
        chargefw.read_parameter_set(write_eem_set(tmp_path / "taken.json", bundled_id, 3.0))

    malformed = tmp_path / "malformed.json"
    malformed.write_text('{"metadata": {"id": "broken"}}')
    with pytest.raises(chargefw.ParseError):
        chargefw.read_parameter_set(malformed)

    unknown_method = write_eem_set(tmp_path / "unknown.json", "water-unknown", 2.0)
    unknown_method.write_text(unknown_method.read_text().replace('"eem"', '"no-such-method"'))
    with pytest.raises(chargefw.ParseError, match="no-such-method"):
        chargefw.read_parameter_set(unknown_method)

    with pytest.raises(FileNotFoundError):
        chargefw.read_parameter_set(tmp_path / "missing.json")


def test_evaluator_with_loaded_parameters_preserves_conformer_assignment_order(
    tmp_path: Any,
) -> None:
    parameter_set = chargefw.read_parameter_set(
        write_eem_set(tmp_path / "set.json", "water-native", 3.0)
    )
    coordinates = water(conformers=2).coordinates.copy()
    coordinates[1] *= 1.5
    molecule = chargefw.Molecule([8, 1, 1], bonds=[[0, 1, 1], [0, 2, 1]], coordinates=coordinates)
    plan = chargefw.assess(
        molecule, method="eem", parameter_set=parameter_set, execution="full"
    ).default_plan
    assert plan is not None
    evaluator = Evaluator(plan, threads=1)

    # The set's own values reproduce a normal calculation, with both conformers concatenated.
    expected = plan.calculate(threads=1)
    assert not np.allclose(expected.assignments[0].values, expected.assignments[1].values)
    np.testing.assert_array_equal(
        evaluator.charges(evaluator.values),
        np.concatenate([assignment.values for assignment in expected.assignments]),
    )
