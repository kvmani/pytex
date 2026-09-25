"""Plans, the unknown mounting rotation, the uncertainty budget and the prose.

Lanes of section 12 of ``docs/architecture/fib_lamella_planning_foundation.md``
covered here:

- *envelope* --- the ``phi``-sweep feasibility fraction is 1 for ``eps`` inside a
  symmetric rectangular envelope, and matches a closed-form arc computation on
  an asymmetric one;
- *cross-check* --- the residual PyTex's TEM navigation solver finds for the
  mounted lamella equals ``eps*``, which ties this module to the validated
  solver; and the vectorized closed-form sweep agrees with the navigation sweep;
- *golden prose* --- ``describe()`` is pinned for a fixed case, and every number
  it states is present in ``to_json_dict()``.
"""

from __future__ import annotations

import json
import math
import re
from pathlib import Path

import numpy as np
import pytest

from pytex.core.conventions import FrameDomain, Handedness
from pytex.core.frames import ReferenceFrame
from pytex.core.lattice import Lattice, Phase
from pytex.core.orientation import Rotation
from pytex.core.symmetry import SymmetrySpec
from pytex.fib import (
    ChamberGeometry,
    FeasibilityClass,
    LamellaSpec,
    MountModel,
    SecondaryPlane,
    UncertaintyInputs,
    lamella_geometry,
    plan_lamella,
    sweep_mounting_rotation,
    target_orbit,
)
from pytex.tem.stage import DoubleTiltStage, EllipticalEnvelope, RectangularEnvelope

GOLDEN = Path(__file__).parent / "golden" / "fib_lamella_plan_describe.txt"


def _phase(name: str, point_group: str, a: float, c: float | None = None) -> Phase:
    frame = ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)
    lattice = Lattice(a, a, c if c is not None else a, 90.0, 90.0, 90.0, crystal_frame=frame)
    return Phase(
        name,
        lattice=lattice,
        symmetry=SymmetrySpec.from_point_group(point_group, reference_frame=frame),
        crystal_frame=frame,
    )


@pytest.fixture(scope="module")
def cubic() -> Phase:
    return _phase("nickel", "m-3m", 3.52)


@pytest.fixture(scope="module")
def tetragonal() -> Phase:
    return _phase("tetragonal", "4/mmm", 3.99, 4.04)


def _rx(degrees: float) -> np.ndarray:
    return Rotation.from_axis_angle([1.0, 0.0, 0.0], math.radians(degrees)).as_matrix()


def _stage(envelope) -> DoubleTiltStage:
    return DoubleTiltStage(envelope=envelope, name="test holder")


def _arc_predicate(eps_deg, sigma, phi_deg, flip, envelope: RectangularEnvelope) -> np.ndarray:
    """Closed-form reachability of the single member (cos e, 0, sigma sin e).

    Mounted front (flip +1) the member sits at w = (-s sigma sin e, c sigma sin e,
    cos e); back it sits at (s sigma sin e, -c sigma sin e, -cos e), which is exactly
    the reverse of the front form at the same phi. With both senses accepted the
    two branches therefore coincide, and the principal tilt is
    alpha = asin(sigma cos phi sin e), tan beta = sigma sin phi tan e.
    """

    e = math.radians(eps_deg)
    phi = np.radians(phi_deg)
    k = sigma
    alpha = np.degrees(np.arcsin(k * np.cos(phi) * math.sin(e)))
    beta = np.degrees(np.arctan(k * np.sin(phi) * math.tan(e)))
    return (
        (alpha >= envelope.alpha_min_deg)
        & (alpha <= envelope.alpha_max_deg)
        & (beta >= envelope.beta_min_deg)
        & (beta <= envelope.beta_max_deg)
    )


class TestEnvelope:
    @pytest.mark.parametrize("eps", [0.0, 5.0, 17.3, 29.5])
    def test_inside_a_symmetric_rectangle_every_rotation_works(self, tetragonal, eps) -> None:
        orbit = target_orbit(tetragonal, (0, 0, 1))
        g = _rx(90.0 - eps)
        geometry = lamella_geometry(g, orbit)
        sweep = sweep_mounting_rotation(
            g, orbit, geometry, mount=MountModel(phi_samples=360, solver="closed_form")
        )
        assert sweep.fraction_front == 1.0 and sweep.fraction_back == 1.0

    def test_an_asymmetric_rectangle_matches_the_closed_form_arcs(self, tetragonal) -> None:
        envelope = RectangularEnvelope(-10.0, 30.0, -20.0, 25.0)
        orbit = target_orbit(tetragonal, (0, 0, 1))
        eps = 20.0
        g = _rx(90.0 - eps)
        geometry = lamella_geometry(g, orbit)
        sigma = geometry.out_of_plane_sign
        mount = MountModel(phi_samples=720, solver="closed_form")
        sweep = sweep_mounting_rotation(g, orbit, geometry, stage=_stage(envelope), mount=mount)
        # Point for point at the sampled rotations...
        np.testing.assert_array_equal(
            sweep.reachable_front, _arc_predicate(eps, sigma, mount.phi_deg, 1, envelope)
        )
        np.testing.assert_array_equal(
            sweep.reachable_back, _arc_predicate(eps, sigma, mount.phi_deg, -1, envelope)
        )
        # ...and as a measure: the analytic arc fraction on a fine grid.
        fine = np.linspace(0.0, 360.0, 200_001)[:-1]
        exact = float(np.mean(_arc_predicate(eps, sigma, fine, 1, envelope)))
        assert sweep.fraction_front == pytest.approx(exact, abs=1.5 / 720)
        assert 0.0 < exact < 1.0

    def test_the_navigation_solver_agrees_with_the_closed_form(self, cubic) -> None:
        g = Rotation.from_axis_angle([0.3, -0.7, 0.4], 1.3).as_matrix()
        orbit = target_orbit(cubic, (1, 1, 1))
        geometry = lamella_geometry(g, orbit)
        stage = _stage(RectangularEnvelope(-12.0, 12.0, -30.0, 30.0))
        closed = sweep_mounting_rotation(
            g, orbit, geometry, stage=stage, mount=MountModel(phi_samples=24, solver="closed_form")
        )
        navigation = sweep_mounting_rotation(
            g, orbit, geometry, stage=stage, mount=MountModel(phi_samples=24, solver="navigation")
        )
        assert navigation.solver == "navigation"
        np.testing.assert_array_equal(closed.reachable_front, navigation.reachable_front)
        np.testing.assert_array_equal(closed.reachable_back, navigation.reachable_back)
        both = closed.reachable_front
        tilt_closed = np.degrees(
            np.arccos(
                np.cos(np.radians(closed.alpha_front_deg[both]))
                * np.cos(np.radians(closed.beta_front_deg[both]))
            )
        )
        tilt_nav = np.degrees(
            np.arccos(
                np.cos(np.radians(navigation.alpha_front_deg[both]))
                * np.cos(np.radians(navigation.beta_front_deg[both]))
            )
        )
        np.testing.assert_allclose(tilt_closed, tilt_nav, atol=1e-6)

    def test_a_coupled_envelope_is_handled_through_its_margin(self, tetragonal) -> None:
        orbit = target_orbit(tetragonal, (0, 0, 1))
        g = _rx(90.0 - 12.0)
        plan = plan_lamella(
            g,
            (0, 0, 1),
            phase=tetragonal,
            stage=_stage(EllipticalEnvelope(30.0, 30.0)),
            mount=MountModel(phi_samples=36, solver="closed_form"),
            orbit=orbit,
            crosscheck=False,
        )
        assert plan.sweep.fraction == 1.0


class TestCrossCheck:
    @pytest.mark.parametrize("seed", [1, 2, 3, 4])
    def test_the_tem_solver_finds_eps_star(self, cubic, seed) -> None:
        rng = np.random.default_rng(seed)
        axis = rng.normal(size=3)
        g = Rotation.from_axis_angle(
            axis / np.linalg.norm(axis), rng.uniform(0, math.pi)
        ).as_matrix()
        plan = plan_lamella(
            g, (1, 1, 2), phase=cubic, mount=MountModel(phi_samples=8, solver="closed_form")
        )
        assert plan.crosscheck_residual_deg == pytest.approx(plan.eps_deg, abs=1e-6)
        assert not any(flag.code == "crosscheck_mismatch" for flag in plan.risk_flags)


class TestFeasibility:
    def _plan(self, tetragonal, eps: float, **kwargs):
        return plan_lamella(
            _rx(90.0 - eps),
            (0, 0, 1),
            phase=tetragonal,
            mount=MountModel(phi_samples=72, solver="closed_form"),
            crosscheck=False,
            **kwargs,
        )

    def test_a_small_residual_is_guaranteed(self, tetragonal) -> None:
        plan = self._plan(tetragonal, 5.0)
        assert plan.feasibility is FeasibilityClass.GUARANTEED

    def test_a_residual_near_the_limit_is_probabilistic_not_guaranteed(self, tetragonal) -> None:
        plan = self._plan(tetragonal, 29.0)
        assert plan.feasibility is FeasibilityClass.PROBABILISTIC
        assert plan.sweep.fraction == 1.0  # every rotation works at the point estimate...
        assert plan.sweep.conservative_fraction < 1.0  # ...but not at the upper bound

    def test_a_residual_beyond_the_holder_is_unreachable(self, tetragonal) -> None:
        plan = self._plan(tetragonal, 60.0)
        assert plan.feasibility is FeasibilityClass.UNREACHABLE
        assert plan.sweep.fraction == 0.0
        assert not plan.is_feasible

    def test_guaranteed_only_at_the_point_estimate_is_flagged(self, tetragonal) -> None:
        # eps + margin = 27 <= 30, but eps + U + margin = 27 + 4.12 > 30.
        plan = self._plan(tetragonal, 22.0)
        assert plan.sweep.guaranteed_at_point and not plan.sweep.guaranteed_at_upper
        assert plan.feasibility is FeasibilityClass.PROBABILISTIC
        assert "point_estimate_only" in {flag.code for flag in plan.risk_flags}

    def test_the_budget_combines_in_quadrature_and_the_verdict_uses_the_upper_bound(
        self, tetragonal
    ) -> None:
        inputs = UncertaintyInputs(
            ebsd_accuracy_deg=0.5, mount_repeatability_deg=2.0, coverage_factor=2.0
        )
        plan = self._plan(tetragonal, 10.0, uncertainty=inputs, grain_spread_deg=1.5)
        expected = math.sqrt(0.5**2 + 1.5**2 + 0.0**2 + 2.0**2)
        assert plan.budget.combined_deg == pytest.approx(expected)
        assert plan.budget.eps_upper_deg == pytest.approx(10.0 + 2.0 * expected)
        assert [item.name for item in plan.budget.components] == [
            "EBSD orientation accuracy",
            "intragranular spread",
            "registration residual",
            "stage and mount repeatability",
        ]

    def test_a_polar_phase_reports_different_flip_branches(self) -> None:
        polar = _phase("polar", "4mm", 3.9, 4.1)
        envelope = RectangularEnvelope(-5.0, 30.0, -30.0, 30.0)
        plan = plan_lamella(
            _rx(90.0 - 20.0),
            (0, 0, 1),
            phase=polar,
            both_senses=False,
            stage=_stage(envelope),
            mount=MountModel(phi_samples=72, solver="navigation"),
        )
        assert plan.sweep.solver == "closed_form"  # one sense only: the closed form is exact
        assert plan.sweep.fraction_front != plan.sweep.fraction_back
        assert "beam_sense_matters" in {flag.code for flag in plan.risk_flags}

    def test_an_uncalibrated_chamber_is_a_critical_risk(self, tetragonal) -> None:
        plan = self._plan(tetragonal, 5.0)
        assert any(
            flag.code == "uncalibrated_azimuth" and flag.severity == "critical"
            for flag in plan.risk_flags
        )
        calibrated = self._plan(
            tetragonal, 5.0, chamber=ChamberGeometry(calibrated=True, calibration_residual_deg=0.3)
        )
        assert "uncalibrated_azimuth" not in {flag.code for flag in calibrated.risk_flags}

    def test_the_subsurface_assumption_is_always_reported(self, tetragonal) -> None:
        plan = self._plan(tetragonal, 5.0)
        assert "subsurface_assumed" in {flag.code for flag in plan.risk_flags}
        assert "not yet against a lamella cut on a real instrument" in plan.describe()

    def test_no_single_alpha_beta_pair_is_reported_as_if_phi_were_known(self, tetragonal) -> None:
        payload = self._plan(tetragonal, 15.0).to_json_dict()
        assert "alpha_deg" not in payload and "beta_deg" not in payload
        assert len(payload["sweep"]["alpha_front_deg"]) == 72


class TestSecondaryScore:
    def test_a_crystal_plane_containing_the_axis_is_edge_on(self, cubic) -> None:
        plan = plan_lamella(
            np.eye(3),
            (1, 0, 0),
            phase=cubic,
            secondary=SecondaryPlane("cube plane", crystal_indices=(0, 0, 1)),
            mount=MountModel(phi_samples=8, solver="closed_form"),
            crosscheck=False,
        )
        assert plan.secondary is not None
        assert plan.secondary.angle_deg == pytest.approx(0.0, abs=1e-9)
        # It does not change the primary answer.
        plain = plan_lamella(
            np.eye(3),
            (1, 0, 0),
            phase=cubic,
            mount=MountModel(phi_samples=8, solver="closed_form"),
            crosscheck=False,
        )
        assert plain.eps_deg == plan.eps_deg
        assert plain.geometry.theta_sample_deg == plan.geometry.theta_sample_deg

    def test_a_sample_plane_is_measured_against_the_beam(self, cubic) -> None:
        plan = plan_lamella(
            np.eye(3),
            (1, 0, 0),
            phase=cubic,
            secondary=SecondaryPlane("boundary", sample_normal=(1.0, 0.0, 0.0)),
            mount=MountModel(phi_samples=8, solver="closed_form"),
            crosscheck=False,
        )
        assert plan.secondary is not None
        assert plan.secondary.angle_deg == pytest.approx(90.0)

    def test_exactly_one_plane_description_is_required(self) -> None:
        with pytest.raises(ValueError):
            SecondaryPlane("both", crystal_indices=(1, 0, 0), sample_normal=(1.0, 0.0, 0.0))


class TestSpec:
    @pytest.mark.parametrize(
        "kwargs",
        [
            {"length_um": 0.0},
            {"width_um": 20.0},
            {"final_thickness_nm": 5000.0},
            {"depth_um": -1.0},
        ],
    )
    def test_impossible_lamellae_are_refused(self, kwargs) -> None:
        with pytest.raises(ValueError):
            LamellaSpec(**kwargs)


def _golden_plan(cubic: Phase):
    g = Rotation.from_bunge_euler(30.0, 40.0, 50.0).as_matrix()
    return plan_lamella(
        g,
        (0, 1, 1),
        phase=cubic,
        grain_id=7,
        location_um=(12.5, 30.0),
        grain_spread_deg=0.4,
        mount=MountModel(phi_samples=36, solver="closed_form"),
    )


class TestProse:
    def test_describe_matches_the_golden_text(self, cubic) -> None:
        text = _golden_plan(cubic).describe()
        # To regenerate after a deliberate wording change, write `text` to GOLDEN
        # and review the diff: the prose is an output, reviewed like any other.
        assert GOLDEN.exists(), f"missing golden file {GOLDEN}"
        assert text == GOLDEN.read_text(encoding="utf-8").rstrip("\n")

    def test_every_number_in_describe_is_in_the_json(self, cubic) -> None:
        plan = _golden_plan(cubic)
        text = plan.describe()
        payload = plan.to_json_dict()
        values: list[float] = []

        def walk(item) -> None:
            if isinstance(item, bool) or item is None:
                return
            if isinstance(item, int | float):
                values.append(float(item))
            elif isinstance(item, dict):
                for value in item.values():
                    walk(value)
            elif isinstance(item, list):
                for value in item:
                    walk(value)
            elif isinstance(item, str):
                for token in re.findall(r"-?\d+(?:\.\d+)?", item):
                    values.append(float(token))

        walk(payload)
        candidates = np.array(values + [100.0 * v for v in values])
        # Numbers inside Miller brackets and citations are text, not measurements.
        stripped = re.sub(r"[<\[(]-?\d[\d -]*[>\])]", " ", text.split(" Sources:")[0])
        stripped = re.sub(r"\d+%", " ", stripped)
        for token in re.findall(r"(?<![\w.])-?\d+(?:\.\d+)?", stripped):
            decimals = len(token.split(".")[1]) if "." in token else 0
            number = float(token)
            assert np.any(
                np.abs(np.round(candidates, decimals) - number) < 10.0 ** (-decimals) / 2 + 1e-12
            ), token

    def test_the_payload_is_json_serializable(self, cubic) -> None:
        json.dumps(_golden_plan(cubic).to_json_dict())

    def test_the_work_order_carries_the_caveat_and_the_uncertainty(self, cubic) -> None:
        rows = dict(_golden_plan(cubic).work_order_lines())
        assert "UNCALIBRATED" in rows["Calibration"]
        assert "+/-" in rows["Expected TEM residual tilt"]
        assert "Validation status" in rows
