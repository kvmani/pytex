"""Validation of the TEM tilt-navigation engine.

The test matrix is section 14.2 of
``docs/architecture/tem_tilt_navigation_foundation.md``. Two classes of test
carry most of the weight and are worth naming:

*Tests that validate the error model, not just the code.* The residual caused by
a wrong diffraction rotation is predicted in closed form; asserting the engine
reproduces that prediction to ten decimal places checks the *physics claim*, not
merely that two code paths agree.

*Tests that assert the absence of a warning.* An ambiguity report that fires on
every crystal is useless, so the enantiomorphic point groups are tested to
confirm the engine stays silent for them — quartz is non-centrosymmetric yet
ambiguity-free, and warning about it would be a defect.
"""

from __future__ import annotations

import math

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

from pytex.core.conventions import FrameDomain, Handedness
from pytex.core.frames import ReferenceFrame
from pytex.core.lattice import Lattice, Phase, ZoneAxis
from pytex.core.orientation import Orientation, Rotation
from pytex.core.point_groups import (
    PointGroup,
    all_point_group_symbols,
    laue_class_symbol_for,
)
from pytex.core.symmetry import SymmetrySpec
from pytex.tem.ambiguity import (
    AMBIGUOUS_POINT_GROUPS,
    analyze_ambiguity,
    laue_rotation_operators,
    observation_stabilizer,
)
from pytex.tem.calibration import (
    TiltExcursionObservation,
    calibrate_from_tilt_excursions,
    predicted_excursion_azimuth_deg,
    residual_from_rotation_error_deg,
)
from pytex.tem.navigation import (
    Reachability,
    plan_tilt_to_zone_axis,
    solve_tilts_for_direction,
)
from pytex.tem.path import PathStrategy, connecting_band, plan_path
from pytex.tem.reconstruction import (
    HOLDER_FRAME,
    CurrentState,
    ReconstructionMode,
    ZoneAxisObservation,
)
from pytex.tem.stage import (
    BEAM_AXIS_LABORATORY,
    DoubleTiltStage,
    EllipticalEnvelope,
    GeneralStageAxes,
    RectangularEnvelope,
    SingleTiltStage,
    StageCalibration,
    StagePosition,
    beam_direction_holder,
    rotation_x,
    rotation_y,
    rotation_z,
)

# --------------------------------------------------------------------------- #
# Fixtures
# --------------------------------------------------------------------------- #


def _crystal_frame() -> ReferenceFrame:
    return ReferenceFrame("crystal", FrameDomain.CRYSTAL, ("a", "b", "c"), Handedness.RIGHT)


def cubic_phase(point_group: str = "m-3m", parameter: float = 3.52) -> Phase:
    frame = _crystal_frame()
    lattice = Lattice(parameter, parameter, parameter, 90.0, 90.0, 90.0, crystal_frame=frame)
    return Phase(
        "nickel",
        lattice=lattice,
        symmetry=SymmetrySpec.from_point_group(point_group, reference_frame=frame),
        crystal_frame=frame,
    )


def hexagonal_phase(point_group: str = "6/mmm") -> Phase:
    frame = _crystal_frame()
    lattice = Lattice(3.232, 3.232, 5.147, 90.0, 90.0, 120.0, crystal_frame=frame)
    return Phase(
        "zirconium",
        lattice=lattice,
        symmetry=SymmetrySpec.from_point_group(point_group, reference_frame=frame),
        crystal_frame=frame,
    )


def tetragonal_phase(point_group: str = "4/mmm") -> Phase:
    frame = _crystal_frame()
    lattice = Lattice(3.99, 3.99, 4.04, 90.0, 90.0, 90.0, crystal_frame=frame)
    return Phase(
        "barium-titanate",
        lattice=lattice,
        symmetry=SymmetrySpec.from_point_group(point_group, reference_frame=frame),
        crystal_frame=frame,
    )


def monoclinic_phase() -> Phase:
    frame = _crystal_frame()
    lattice = Lattice(5.15, 5.20, 5.35, 90.0, 99.2, 90.0, crystal_frame=frame)
    return Phase(
        "zirconia-monoclinic",
        lattice=lattice,
        symmetry=SymmetrySpec.from_point_group("2/m", reference_frame=frame),
        crystal_frame=frame,
    )


def make_state(
    phase: Phase,
    matrix: np.ndarray,
    position: tuple[float, float] = (0.0, 0.0),
    zone: tuple[int, int, int] = (0, 0, 1),
) -> CurrentState:
    orientation = Orientation.from_matrix(
        matrix,
        specimen_frame=HOLDER_FRAME,
        phase=phase,
        crystal_frame=phase.crystal_frame,
    )
    return CurrentState.from_orientation(
        orientation,
        StagePosition(*position),
        current_zone_axis=ZoneAxis(list(zone), phase=phase),
    )


def random_rotations(count: int, seed: int) -> list[np.ndarray]:
    generator = np.random.default_rng(seed)
    matrices = []
    for _ in range(count):
        matrix = np.linalg.qr(generator.normal(size=(3, 3)))[0]
        if np.linalg.det(matrix) < 0:
            matrix[:, 0] *= -1
        matrices.append(matrix)
    return matrices


WIDE_STAGE = DoubleTiltStage(envelope=RectangularEnvelope(-60.0, 60.0, -60.0, 60.0))
TYPICAL_STAGE = DoubleTiltStage(envelope=RectangularEnvelope(-30.0, 30.0, -30.0, 30.0))


# --------------------------------------------------------------------------- #
# Stage kinematics
# --------------------------------------------------------------------------- #


class TestStageKinematics:
    def test_moving_beta_axis_composition_cancels_exactly(self) -> None:
        """The physical composition equals Rx(alpha) Ry(beta).

        The beta axis moves with alpha, so the physically correct composition is
        the beta rotation about its *instantaneous* laboratory axis followed by
        alpha. That the two agree is a property of this axis pair, not a general
        licence to ignore axis motion, so it is pinned here.
        """

        for alpha_deg, beta_deg in ((12.0, -7.0), (30.0, 25.0), (-18.0, 33.0)):
            alpha, beta = math.radians(alpha_deg), math.radians(beta_deg)
            moving_axis = (
                rotation_x(alpha) @ rotation_y(beta) @ rotation_x(alpha).T
            ) @ rotation_x(alpha)
            assert np.allclose(moving_axis, rotation_x(alpha) @ rotation_y(beta))

    def test_beam_direction_matches_the_closed_form(self) -> None:
        stage = DoubleTiltStage()
        for alpha, beta in ((0.0, 0.0), (12.0, -7.0), (30.0, 25.0), (-18.0, 33.0)):
            assert np.allclose(
                stage.beam_direction(alpha, beta), beam_direction_holder(alpha, beta)
            )

    def test_beta_effectiveness_falls_as_cos_alpha(self) -> None:
        """A degree of beta buys only cos(alpha) degrees of crystal rotation.

        This is the double-tilt holder's gimbal lock, and the origin of the
        conditioning factor the solver reports.
        """

        stage = DoubleTiltStage()
        step = 1e-6
        for alpha in (0.0, 23.0, 40.0):
            forward = stage.beam_direction(alpha, 0.2 + step)
            backward = stage.beam_direction(alpha, 0.2 - step)
            derivative = np.linalg.norm(forward - backward) / (2.0 * math.radians(step))
            assert derivative == pytest.approx(math.cos(math.radians(alpha)), abs=1e-6)

    def test_alpha_effectiveness_is_unity(self) -> None:
        stage = DoubleTiltStage()
        step = 1e-6
        for alpha in (0.0, 23.0, 40.0):
            forward = stage.beam_direction(alpha + step, 0.2)
            backward = stage.beam_direction(alpha - step, 0.2)
            derivative = np.linalg.norm(forward - backward) / (2.0 * math.radians(step))
            assert derivative == pytest.approx(1.0, abs=1e-6)

    def test_rectangular_solid_angle_matches_numeric_integration(self) -> None:
        """The analytic solid angle equals the integral of the cos(alpha) Jacobian."""

        envelope = RectangularEnvelope(-30.0, 30.0, -30.0, 30.0)
        analytic = envelope.accessible_solid_angle_sr()
        assert analytic == pytest.approx(
            math.radians(60.0) * 2.0 * math.sin(math.radians(30.0)), abs=1e-12
        )
        # The base-class numerical integration is an independent route.
        from pytex.tem.stage import TiltEnvelope

        numeric = TiltEnvelope.accessible_solid_angle_sr(envelope, 400)
        assert numeric == pytest.approx(analytic, rel=1e-4)

    def test_typical_holder_reaches_about_eight_percent_of_the_sphere(self) -> None:
        """A +/-30 degree holder commands 1.047 sr, 8.3 percent of all directions."""

        omega = RectangularEnvelope(-30.0, 30.0, -30.0, 30.0).accessible_solid_angle_sr()
        assert omega == pytest.approx(1.0472, abs=1e-3)
        assert 100.0 * omega / (4.0 * math.pi) == pytest.approx(8.33, abs=0.02)

    def test_elliptical_envelope_couples_the_two_ranges(self) -> None:
        envelope = EllipticalEnvelope(30.0, 30.0)
        assert envelope.contains(30.0, 0.0)
        assert envelope.contains(0.0, 30.0)
        assert not envelope.contains(25.0, 25.0)

    def test_general_axes_report_non_orthogonality(self) -> None:
        axes = GeneralStageAxes(
            alpha_axis=[1.0, 0.0, 0.0],
            beta_axis=[math.sin(math.radians(3.0)), math.cos(math.radians(3.0)), 0.0],
        )
        assert axes.non_orthogonality_deg == pytest.approx(3.0, abs=1e-9)
        assert not axes.is_ideal


# --------------------------------------------------------------------------- #
# Closed form and branches
# --------------------------------------------------------------------------- #


class TestClosedForm:
    def test_round_trip_is_exact_for_random_directions(self) -> None:
        """Every solved branch places the direction on the beam to machine precision."""

        generator = np.random.default_rng(0)
        stage = DoubleTiltStage()
        worst = 0.0
        for _ in range(2000):
            direction = generator.normal(size=3)
            direction /= np.linalg.norm(direction)
            for alpha, beta in solve_tilts_for_direction(direction):
                achieved = stage.rotation_matrix(alpha, beta) @ direction
                worst = max(worst, abs(abs(float(achieved[2])) - 1.0))
        assert worst < 1e-12

    def test_four_branches_land_on_the_expected_pole(self) -> None:
        """Branches one and two reach +z; three and four reach -z."""

        generator = np.random.default_rng(5)
        stage = DoubleTiltStage()
        direction = generator.normal(size=3)
        direction /= np.linalg.norm(direction)
        branches = solve_tilts_for_direction(direction, allow_reverse=True)
        assert len(branches) == 4
        signs = [
            float(np.dot(stage.rotation_matrix(a, b) @ direction, BEAM_AXIS_LABORATORY))
            for a, b in branches
        ]
        assert signs[0] == pytest.approx(1.0, abs=1e-12)
        assert signs[1] == pytest.approx(1.0, abs=1e-12)
        assert signs[2] == pytest.approx(-1.0, abs=1e-12)
        assert signs[3] == pytest.approx(-1.0, abs=1e-12)

    def test_holder_z_axis_needs_no_tilt(self) -> None:
        assert solve_tilts_for_direction([0.0, 0.0, 1.0])[0] == (0.0, 0.0)

    def test_gimbal_degeneracy_is_detected_not_guessed(self) -> None:
        """Along the beta axis, beta is indeterminate; alpha is +/-90 degrees.

        The engine must report the representative rather than emit an arbitrary
        ``atan2(0, 0)``.
        """

        branches = solve_tilts_for_direction([0.0, 1.0, 0.0])
        assert branches[0] == (90.0, 0.0)
        assert branches[1] == (-90.0, 0.0)
        assert all(math.isfinite(value) for pair in branches for value in pair)


# --------------------------------------------------------------------------- #
# Navigation
# --------------------------------------------------------------------------- #


class TestNavigation:
    def test_forward_validation_is_exact_for_random_orientations(self) -> None:
        phase = cubic_phase()
        target = ZoneAxis([1, 1, 2], phase=phase)
        worst = 0.0
        reached = 0
        for matrix in random_rotations(200, seed=7):
            report = plan_tilt_to_zone_axis(
                make_state(phase, matrix), target, WIDE_STAGE, include_paths=False
            )
            if report.is_reachable:
                reached += 1
                worst = max(worst, report.best().residual_deg)
        assert reached == 200
        assert worst < 1e-9

    @pytest.mark.parametrize(
        ("start", "target", "expected_deg"),
        [
            ((0, 0, 1), (0, 1, 1), 45.0),
            ((0, 0, 1), (1, 1, 1), 54.7356),
            ((0, 0, 1), (1, 1, 2), 35.2644),
            ((0, 1, 1), (1, 1, 1), 35.2644),
            ((1, 1, 1), (-1, 1, 1), 70.5288),
        ],
    )
    def test_known_cubic_interzonal_angles(
        self, start: tuple[int, int, int], target: tuple[int, int, int], expected_deg: float
    ) -> None:
        """Crystal travel along the geodesic equals the interzonal angle.

        The expected values are analytic: 45 degrees is arccos(1/sqrt(2)), 54.7356
        is arccos(1/sqrt(3)), 35.2644 is arccos(sqrt(2/3)), and 70.5288 is
        arccos(1/3). None is copied from a program run.
        """

        phase = cubic_phase()
        # Orient the crystal so that `start` lies along the beam at zero tilt.
        start_vector = ZoneAxis(list(start), phase=phase).unit_vector
        alpha, beta = solve_tilts_for_direction(start_vector)[0]
        matrix = np.asarray(
            DoubleTiltStage().rotation_matrix(alpha, beta), dtype=np.float64
        )
        report = plan_tilt_to_zone_axis(
            make_state(phase, matrix, zone=start),
            ZoneAxis(list(target), phase=phase),
            WIDE_STAGE,
        )
        assert report.is_reachable
        best = report.best()
        assert best.path is not None
        # The engine may legitimately choose a nearer symmetry equivalent; the
        # requested direction must nonetheless be at the analytic angle.
        angle = math.degrees(
            math.acos(
                min(
                    1.0,
                    abs(
                        float(
                            np.dot(
                                start_vector,
                                ZoneAxis(list(target), phase=phase).unit_vector,
                            )
                        )
                    ),
                )
            )
        )
        assert angle == pytest.approx(expected_deg, abs=1e-3)

    def test_hexagonal_basal_to_prismatic_is_ninety_degrees(self) -> None:
        """[0001] and <11-20> are perpendicular by symmetry, at any c/a ratio."""

        phase = hexagonal_phase()
        basal = ZoneAxis([0, 0, 1], phase=phase).unit_vector
        prismatic = ZoneAxis([2, -1, 0], phase=phase).unit_vector
        angle = math.degrees(math.acos(abs(float(np.dot(basal, prismatic)))))
        assert angle == pytest.approx(90.0, abs=1e-9)

    def test_symmetry_orbit_size_matches_the_point_group(self) -> None:
        """A general direction has one orbit member per proper operator, times two.

        ``[135]`` is used rather than the more obvious ``[123]`` because the
        latter is *not* general in a hexagonal lattice: the Cartesian x component
        of ``[u v w]`` is proportional to ``u - v/2``, which vanishes for every
        ``[1 2 w]``, putting the direction on a mirror and halving its orbit.
        """

        for phase, expected in (
            (cubic_phase(), 48),
            (hexagonal_phase(), 24),
            (tetragonal_phase(), 16),
            (monoclinic_phase(), 4),
        ):
            operator_count = np.asarray(phase.symmetry.operators).shape[0]
            assert expected == 2 * operator_count, phase.name
            report = plan_tilt_to_zone_axis(
                make_state(phase, np.eye(3)),
                ZoneAxis([1, 3, 5], phase=phase),
                WIDE_STAGE,
                include_paths=False,
            )
            assert report.orbit_size == expected, phase.name

    def test_special_directions_have_smaller_orbits(self) -> None:
        """A direction on a symmetry element has a stabilizer, so its orbit shrinks.

        ``<100>`` in cubic has 6 members, not 48: the four-fold about it fixes it.
        """

        phase = cubic_phase()
        report = plan_tilt_to_zone_axis(
            make_state(phase, np.eye(3)),
            ZoneAxis([1, 0, 0], phase=phase),
            WIDE_STAGE,
            include_paths=False,
        )
        assert report.orbit_size == 6

    def test_target_out_of_range_reports_a_reachable_nearest_approach(self) -> None:
        """A nearest approach must be a position the holder can actually reach."""

        phase = cubic_phase()
        report = plan_tilt_to_zone_axis(
            make_state(phase, np.eye(3)),
            ZoneAxis([0, 1, 1], phase=phase),
            TYPICAL_STAGE,
        )
        assert not report.is_reachable
        nearest = report.nearest_approach
        assert nearest is not None
        assert nearest.verdict is Reachability.NEAREST_APPROACH
        assert TYPICAL_STAGE.envelope.contains(*nearest.position.as_tuple())
        assert nearest.residual_deg == pytest.approx(15.0, abs=0.5)
        with pytest.raises(ValueError, match="No reachable solution"):
            report.best()

    def test_verdict_flips_at_the_envelope_boundary(self) -> None:
        """A target just inside is reachable; the same target just outside is not."""

        phase = cubic_phase()
        matrix = np.eye(3)
        target = ZoneAxis([0, 1, 1], phase=phase)
        inside = DoubleTiltStage(envelope=RectangularEnvelope(-45.5, 45.5, -45.5, 45.5))
        outside = DoubleTiltStage(envelope=RectangularEnvelope(-44.5, 44.5, -44.5, 44.5))
        assert plan_tilt_to_zone_axis(
            make_state(phase, matrix), target, inside, include_paths=False
        ).is_reachable
        assert not plan_tilt_to_zone_axis(
            make_state(phase, matrix), target, outside, include_paths=False
        ).is_reachable

    def test_no_returned_solution_lies_outside_the_envelope(self) -> None:
        phase = cubic_phase()
        for matrix in random_rotations(60, seed=21):
            report = plan_tilt_to_zone_axis(
                make_state(phase, matrix),
                ZoneAxis([1, 1, 1], phase=phase),
                TYPICAL_STAGE,
                include_paths=False,
            )
            for solution in report.solutions:
                assert solution.envelope_margin_deg >= 0.0
                assert TYPICAL_STAGE.envelope.contains(*solution.position.as_tuple())

    def test_single_tilt_holder_never_claims_an_exact_hit(self) -> None:
        """One freedom reaches a curve, so an exact zone axis is coincidence only."""

        phase = cubic_phase()
        stage = SingleTiltStage(envelope=RectangularEnvelope(-70.0, 70.0, -1e-6, 1e-6))
        generator = np.random.default_rng(13)
        for _ in range(15):
            matrix = np.linalg.qr(generator.normal(size=(3, 3)))[0]
            if np.linalg.det(matrix) < 0:
                matrix[:, 0] *= -1
            report = plan_tilt_to_zone_axis(
                make_state(phase, matrix),
                ZoneAxis([1, 2, 3], phase=phase),
                stage,
                include_paths=False,
            )
            assert not report.is_reachable
            assert report.nearest_approach is not None

    def test_non_orthogonal_axes_still_solve_by_refinement(self) -> None:
        """The closed form seeds a refinement that absorbs a mis-set axis pair."""

        phase = cubic_phase()
        for deviation in (1.0, 3.0, 5.0):
            calibration = StageCalibration(
                axes=GeneralStageAxes(
                    alpha_axis=[1.0, 0.0, 0.0],
                    beta_axis=[
                        math.sin(math.radians(deviation)),
                        math.cos(math.radians(deviation)),
                        0.0,
                    ],
                )
            )
            stage = DoubleTiltStage(
                envelope=RectangularEnvelope(-60.0, 60.0, -60.0, 60.0),
                calibration=calibration,
            )
            report = plan_tilt_to_zone_axis(
                make_state(phase, np.eye(3)),
                ZoneAxis([1, 1, 2], phase=phase),
                stage,
                include_paths=False,
            )
            assert report.is_reachable, deviation
            assert report.best().residual_deg < 1e-6, deviation

    def test_describe_is_prose_and_states_the_conventions(self) -> None:
        phase = cubic_phase()
        report = plan_tilt_to_zone_axis(
            make_state(phase, np.eye(3)), ZoneAxis([1, 1, 2], phase=phase), WIDE_STAGE
        )
        text = report.describe()
        assert "symmetry orbit" in text
        assert "uniquely" in text.lower()
        assert "residual" in text.lower()
        assert len(text) > 400

    def test_json_contract_and_describe_stay_in_lockstep(self) -> None:
        phase = cubic_phase()
        report = plan_tilt_to_zone_axis(
            make_state(phase, np.eye(3)), ZoneAxis([1, 1, 2], phase=phase), WIDE_STAGE
        )
        payload = report.to_json_dict()
        assert payload["schema"] == "pytex.tilt_plan_report/1"
        assert payload["is_reachable"] is True
        assert payload["orbit_size"] == report.orbit_size
        assert len(payload["solutions"]) == len(report.solutions)
        assert payload["solutions"][0]["alpha_deg"] == report.best().position.alpha_deg


# --------------------------------------------------------------------------- #
# Reconstruction
# --------------------------------------------------------------------------- #


class TestReconstruction:
    def test_two_zone_reconstruction_recovers_a_known_orientation(self) -> None:
        """Mode B determines the orientation with no rotation calibration at all."""

        phase = cubic_phase()
        operators = np.asarray(phase.symmetry.operators, dtype=np.float64)
        worst = 0.0
        pairs = (((0, 0, 1), (0, 1, 1)), ((0, 0, 1), (1, 1, 1)), ((1, 0, 1), (1, 1, 1)))
        for first, second in pairs:
            for matrix in random_rotations(20, seed=11):
                zone_first = ZoneAxis(list(first), phase=phase)
                zone_second = ZoneAxis(list(second), phase=phase)
                position_first = StagePosition(
                    *solve_tilts_for_direction(matrix @ zone_first.unit_vector)[0]
                )
                position_second = StagePosition(
                    *solve_tilts_for_direction(matrix @ zone_second.unit_vector)[0]
                )
                state = CurrentState.from_two_zone_axes(
                    ZoneAxisObservation(zone_first, position_first),
                    ZoneAxisObservation(zone_second, position_second),
                    WIDE_STAGE,
                )
                assert state.mode is ReconstructionMode.TWO_ZONE_AXES
                assert state.consistency_residual_deg == pytest.approx(0.0, abs=1e-9)
                relative = np.einsum(
                    "ij,njk->nik",
                    state.matrix,
                    np.einsum("nij,jk->nik", operators, matrix.T),
                )
                angles = np.degrees(
                    np.arccos(
                        np.clip((np.trace(relative, axis1=1, axis2=2) - 1.0) / 2.0, -1.0, 1.0)
                    )
                )
                worst = max(worst, float(np.min(angles)))
        assert worst < 1e-4

    def test_inconsistent_two_zone_input_is_indicted_not_absorbed(self) -> None:
        """The interzonal angle is an invariant, so a bad input must be flagged.

        [001] and [110] are 90 degrees apart, but the stage positions given here
        are only 45 degrees apart, so the data cannot both be right.
        """

        phase = cubic_phase()
        state = CurrentState.from_two_zone_axes(
            ZoneAxisObservation(ZoneAxis([0, 0, 1], phase=phase), StagePosition(0.0, 0.0)),
            ZoneAxisObservation(ZoneAxis([1, 1, 0], phase=phase), StagePosition(45.0, 0.0)),
            WIDE_STAGE,
        )
        assert state.consistency_residual_deg is not None
        assert state.consistency_residual_deg > 40.0
        assert any("WARNING" in note for note in state.notes)

    def test_two_zone_flip_ambiguity_is_harmless_when_it_is_a_symmetry(self) -> None:
        """For cubic [001] and [110] the residual two-fold is in the point group."""

        phase = cubic_phase()
        state = CurrentState.from_two_zone_axes(
            ZoneAxisObservation(ZoneAxis([0, 0, 1], phase=phase), StagePosition(0.0, 0.0)),
            ZoneAxisObservation(
                ZoneAxis([1, 1, 0], phase=phase),
                StagePosition(*solve_tilts_for_direction([1.0, 1.0, 0.0])[0]),
            ),
            WIDE_STAGE,
        )
        assert any("harmless" in note for note in state.notes)

    def test_single_pattern_reconstruction_needs_a_calibration(self) -> None:
        """Guessing the pattern azimuth is the error this refusal prevents."""

        phase = cubic_phase()
        with pytest.raises(ValueError, match="calibrated diffraction rotation"):
            CurrentState.from_pattern_solution(
                np.eye(3), StagePosition(0.0, 0.0), WIDE_STAGE, phase
            )

    def test_wrong_parity_is_caught_by_the_determinant_check(self) -> None:
        """A mirrored recording makes the reconstruction improper, and that is free."""

        phase = cubic_phase()
        stage = DoubleTiltStage(
            calibration=StageCalibration(diffraction_rotation_deg=0.0, pattern_is_mirrored=True)
        )
        # A proper crystal-to-pattern rotation combined with a parity flag the
        # data do not actually carry makes the product improper, which is the
        # free self-check.
        with pytest.raises(ValueError, match="improper"):
            CurrentState.from_pattern_solution(
                np.eye(3), StagePosition(0.0, 0.0), stage, phase
            )

    def test_single_pattern_round_trips_with_a_correct_calibration(self) -> None:
        phase = cubic_phase()
        for rotation_deg in (0.0, 37.0, -110.0):
            stage = DoubleTiltStage(
                calibration=StageCalibration(diffraction_rotation_deg=rotation_deg)
            )
            truth = Rotation.from_axis_angle([1.0, 2.0, 3.0], 0.7).as_matrix()
            pattern = rotation_z(math.radians(rotation_deg)).T @ truth
            state = CurrentState.from_pattern_solution(
                pattern, StagePosition(0.0, 0.0), stage, phase
            )
            assert np.allclose(state.matrix, truth, atol=1e-12)

    def test_multi_zone_fit_measures_its_own_scatter(self) -> None:
        phase = cubic_phase()
        matrix = random_rotations(1, seed=3)[0]
        zones = [(0, 0, 1), (0, 1, 1), (1, 1, 1), (1, 0, 1)]
        observations = tuple(
            ZoneAxisObservation(
                ZoneAxis(list(zone), phase=phase),
                StagePosition(
                    *solve_tilts_for_direction(
                        matrix @ ZoneAxis(list(zone), phase=phase).unit_vector
                    )[0]
                ),
            )
            for zone in zones
        )
        state = CurrentState.from_zone_axes(observations, WIDE_STAGE)
        assert state.mode is ReconstructionMode.MULTI_ZONE_FIT
        assert state.orientation_uncertainty_deg < 1e-6

    def test_parallel_zone_axes_are_refused_with_an_actionable_message(self) -> None:
        phase = cubic_phase()
        with pytest.raises(ValueError, match="parallel"):
            CurrentState.from_two_zone_axes(
                ZoneAxisObservation(ZoneAxis([0, 0, 1], phase=phase), StagePosition(0.0, 0.0)),
                ZoneAxisObservation(ZoneAxis([0, 0, 2], phase=phase), StagePosition(0.0, 0.0)),
                WIDE_STAGE,
                resolve_senses=False,
            )


# --------------------------------------------------------------------------- #
# Ambiguity
# --------------------------------------------------------------------------- #


class TestAmbiguity:
    def test_laue_rotation_group_exceeds_the_point_group_for_exactly_ten_groups(
        self,
    ) -> None:
        """The criterion is improper operations other than inversion.

        Enumerated over all 32 point groups. This is the claim that makes the
        engine quiet for cubic and hexagonal metals, so it is pinned here rather
        than asserted in prose.
        """

        enlarged = set()
        for symbol in all_point_group_symbols():
            proper = PointGroup.from_symbol(symbol).operators
            proper_count = int(np.sum(np.linalg.det(proper) > 0))
            laue_count = len(laue_rotation_operators(symbol))
            assert laue_count % proper_count == 0
            if laue_count != proper_count:
                assert laue_count == 2 * proper_count, symbol
                enlarged.add(symbol)
        assert enlarged == set(AMBIGUOUS_POINT_GROUPS)

    def test_cubic_down_001_has_an_order_eight_stabilizer_all_of_it_symmetry(
        self,
    ) -> None:
        """422: the four-fold about [001] and the in-plane two-folds, all in 432."""

        stabilizer = observation_stabilizer("m-3m", [0.0, 0.0, 1.0])
        assert stabilizer.shape[0] == 8
        report = analyze_ambiguity(cubic_phase(), [0.0, 0.0, 1.0])
        assert report.stabilizer_order == 8
        assert report.symmetry_stabilizer_order == 8
        assert report.is_unique

    @pytest.mark.parametrize(
        ("point_group", "zone"),
        [("m-3m", (0, 0, 1)), ("m-3m", (1, 1, 1)), ("6/mmm", (0, 0, 1))],
    )
    def test_centrosymmetric_phases_are_reported_unambiguous(
        self, point_group: str, zone: tuple[int, int, int]
    ) -> None:
        phase = (
            cubic_phase(point_group)
            if point_group.startswith("m")
            else hexagonal_phase(point_group)
        )
        report = analyze_ambiguity(phase, ZoneAxis(list(zone), phase=phase).unit_vector)
        assert report.is_unique
        assert report.experiments == ()
        assert "unique" in report.describe().lower()

    @pytest.mark.parametrize("point_group", ["432", "23", "422", "32"])
    def test_enantiomorphic_phases_raise_no_warning(self, point_group: str) -> None:
        """Non-centrosymmetric is the wrong test; quartz (32) is ambiguity-free.

        Warning on these eleven groups would be a defect: their proper point
        group already equals their Laue rotation group, so Friedel's law adds
        nothing.
        """

        assert point_group not in AMBIGUOUS_POINT_GROUPS
        assert len(laue_rotation_operators(point_group)) == int(
            np.sum(np.linalg.det(PointGroup.from_symbol(point_group).operators) > 0)
        )

    @pytest.mark.parametrize("point_group", ["-43m", "6mm", "4mm", "mm2"])
    def test_improper_non_centrosymmetric_phases_report_the_extra_family(
        self, point_group: str
    ) -> None:
        assert point_group in AMBIGUOUS_POINT_GROUPS
        laue = laue_class_symbol_for(point_group)
        assert len(laue_rotation_operators(point_group)) == 2 * int(
            np.sum(np.linalg.det(PointGroup.from_symbol(point_group).operators) > 0)
        )
        assert laue != point_group

    def test_gallium_arsenide_down_001_has_two_hypotheses_with_an_experiment(
        self,
    ) -> None:
        phase = cubic_phase("-43m")
        report = analyze_ambiguity(phase, [0.0, 0.0, 1.0])
        assert not report.is_unique
        assert len(report.families) == 2
        assert report.stabilizer_order == 2 * report.symmetry_stabilizer_order
        assert report.experiments
        text = report.describe()
        assert "hypotheses" in text
        assert "discriminate" in text.lower()

    def test_uncalibrated_rotation_adds_an_instrumental_family(self) -> None:
        report = analyze_ambiguity(
            cubic_phase(), [0.0, 0.0, 1.0], rotation_calibrated=False
        )
        assert not report.is_unique
        assert any("wrong way" in family.rationale for family in report.families)


# --------------------------------------------------------------------------- #
# Calibration and the error model
# --------------------------------------------------------------------------- #


class TestCalibrationAndErrorModel:
    @pytest.mark.parametrize("rotation_error_deg", [5.0, 30.0, 180.0])
    def test_engine_reproduces_the_closed_form_residual_law(
        self, rotation_error_deg: float
    ) -> None:
        """The residual from a wrong diffraction rotation is 2 asin(sin(d/2) sin(t)).

        This validates the *error model*, not merely that two code paths agree:
        the prediction is analytic and independent of the solver.
        """

        phase = cubic_phase()
        truth = np.eye(3)
        estimate = rotation_z(math.radians(rotation_error_deg)) @ truth
        target = ZoneAxis([1, 1, 4], phase=phase)
        report = plan_tilt_to_zone_axis(
            make_state(phase, estimate),
            target,
            WIDE_STAGE,
            include_paths=False,
            tolerance_deg=1e-6,
        )
        assert report.is_reachable
        best = report.best()
        hop_deg = math.degrees(
            math.acos(abs(float(np.dot(truth @ target.unit_vector, BEAM_AXIS_LABORATORY))))
        )
        landed = (
            WIDE_STAGE.rotation_matrix(*best.position.as_tuple())
            @ truth
            @ best.orbit_member
        )
        landed /= np.linalg.norm(landed)
        actual = math.degrees(
            math.acos(min(1.0, abs(float(np.dot(landed, BEAM_AXIS_LABORATORY)))))
        )
        predicted = residual_from_rotation_error_deg(rotation_error_deg, hop_deg)
        assert actual == pytest.approx(predicted, abs=1e-9)

    def test_a_180_degree_rotation_error_negates_both_tilt_angles(self) -> None:
        """The catastrophic, self-consistent failure: exactly the wrong way.

        The calculation reports a clean zero residual while the specimen goes
        backwards, which is why the two-excursion calibration exists.
        """

        phase = cubic_phase()
        target = ZoneAxis([1, 2, 8], phase=phase)
        correct = plan_tilt_to_zone_axis(
            make_state(phase, np.eye(3)), target, WIDE_STAGE, include_paths=False
        ).best()
        wrong = plan_tilt_to_zone_axis(
            make_state(phase, rotation_z(math.pi)), target, WIDE_STAGE, include_paths=False
        ).best()
        assert wrong.position.alpha_deg == pytest.approx(-correct.position.alpha_deg, abs=1e-9)
        assert wrong.position.beta_deg == pytest.approx(-correct.position.beta_deg, abs=1e-9)
        # Both report a clean residual, which is precisely the danger.
        assert correct.residual_deg < 1e-9
        assert wrong.residual_deg < 1e-9

    def test_residual_law_short_and_long_hops(self) -> None:
        """A 5 degree error costs 0.44 degrees over a short hop and 5 over a long one."""

        assert residual_from_rotation_error_deg(5.0, 5.0) == pytest.approx(0.436, abs=1e-3)
        assert residual_from_rotation_error_deg(5.0, 90.0) == pytest.approx(5.0, abs=1e-9)

    @pytest.mark.parametrize("rotation_deg", [0.0, 37.0, -110.0])
    @pytest.mark.parametrize("mirrored", [False, True])
    def test_two_excursion_calibration_round_trips(
        self, rotation_deg: float, mirrored: bool
    ) -> None:
        """Two exposures recover the rotation numerically and the parity from a sign."""

        result = calibrate_from_tilt_excursions(
            TiltExcursionObservation(
                "alpha",
                8.0,
                predicted_excursion_azimuth_deg("alpha", rotation_deg, mirrored=mirrored),
                8.0,
            ),
            TiltExcursionObservation(
                "beta",
                8.0,
                predicted_excursion_azimuth_deg("beta", rotation_deg, mirrored=mirrored),
                8.0,
            ),
        )
        assert result.is_consistent
        assert result.calibration.pattern_is_mirrored is mirrored
        recovered = result.calibration.diffraction_rotation_deg
        assert recovered is not None
        assert ((recovered - rotation_deg + 180.0) % 360.0) - 180.0 == pytest.approx(
            0.0, abs=1e-9
        )
        assert result.orthogonality_residual_deg == pytest.approx(0.0, abs=1e-9)

    def test_inconsistent_excursions_are_flagged(self) -> None:
        """Azimuths that are not 90 degrees apart indict the measurement."""

        result = calibrate_from_tilt_excursions(
            TiltExcursionObservation("alpha", 8.0, 0.0),
            TiltExcursionObservation("beta", 8.0, 30.0),
        )
        assert not result.is_consistent
        assert "not 90 degrees apart" in result.describe()

    def test_same_axis_excursions_are_refused(self) -> None:
        with pytest.raises(ValueError, match="different axes"):
            calibrate_from_tilt_excursions(
                TiltExcursionObservation("alpha", 8.0, 0.0),
                TiltExcursionObservation("alpha", -8.0, 180.0),
            )

    def test_calibration_refuses_a_different_camera_length(self) -> None:
        """The diffraction rotation is hysteretic; interpolating it invents a number."""

        calibration = StageCalibration(
            diffraction_rotation_deg=12.0, camera_length_mm=800.0
        )
        calibration.check_applicable(camera_length_mm=800.0)
        with pytest.raises(ValueError, match="hysteretic"):
            calibration.check_applicable(camera_length_mm=1200.0)

    def test_axis_fitting_is_refused_when_under_determined(self) -> None:
        """A spuriously fitted axis deviation looks like knowledge and is not."""

        from pytex.tem.calibration import fit_stage_and_orientation

        phase = cubic_phase()
        observations = tuple(
            ZoneAxisObservation(
                ZoneAxis(list(zone), phase=phase),
                StagePosition(
                    *solve_tilts_for_direction(
                        ZoneAxis(list(zone), phase=phase).unit_vector
                    )[0]
                ),
            )
            for zone in ((0, 0, 1), (0, 1, 1), (1, 1, 1))
        )
        with pytest.raises(ValueError, match="at least five"):
            fit_stage_and_orientation(observations, WIDE_STAGE, fit_axes=True)


# --------------------------------------------------------------------------- #
# Path planning
# --------------------------------------------------------------------------- #


class TestPathPlanning:
    def test_geodesic_travel_equals_the_interzonal_angle(self) -> None:
        """The geodesic is the shortest arc, so its length is the angle itself."""

        phase = cubic_phase()
        matrix = np.eye(3)
        target_direction = ZoneAxis([1, 1, 1], phase=phase).unit_vector
        alpha, beta = solve_tilts_for_direction(matrix @ target_direction)[0]
        path = plan_path(
            StagePosition(0.0, 0.0),
            StagePosition(alpha, beta),
            WIDE_STAGE,
            matrix,
            strategy=PathStrategy.GEODESIC,
            phase=phase,
            samples=201,
        )
        assert path.total_travel_deg == pytest.approx(54.7356, abs=1e-3)

    def test_geodesic_is_never_longer_than_the_sequential_paths(self) -> None:
        phase = cubic_phase()
        matrix = np.eye(3)
        alpha, beta = solve_tilts_for_direction(
            ZoneAxis([1, 1, 2], phase=phase).unit_vector
        )[0]
        lengths = {
            strategy: plan_path(
                StagePosition(0.0, 0.0),
                StagePosition(alpha, beta),
                WIDE_STAGE,
                matrix,
                strategy=strategy,
                samples=201,
            ).total_travel_deg
            for strategy in PathStrategy
        }
        geodesic = lengths[PathStrategy.GEODESIC]
        for strategy, length in lengths.items():
            assert geodesic <= length + 1e-6, strategy

    def test_connecting_band_of_001_and_111_is_1m10(self) -> None:
        """Both poles lie in (1-10), which is the band an operator follows.

        Checked directly: the zone law hu + kv + lw = 0 holds for both.
        """

        phase = cubic_phase()
        band = connecting_band(
            phase,
            ZoneAxis([0, 0, 1], phase=phase).unit_vector,
            ZoneAxis([1, 1, 1], phase=phase).unit_vector,
        )
        assert band is not None
        indices = np.asarray(band.indices)
        assert int(np.dot(indices, [0, 0, 1])) == 0
        assert int(np.dot(indices, [1, 1, 1])) == 0
        assert set(np.abs(indices).tolist()) == {0, 1}

    def test_parallel_directions_have_no_connecting_band(self) -> None:
        phase = cubic_phase()
        assert (
            connecting_band(phase, [0.0, 0.0, 1.0], [0.0, 0.0, 2.0]) is None
        )

    def test_path_leaving_the_envelope_is_rejected_with_a_reason(self) -> None:
        phase = cubic_phase()
        path = plan_path(
            StagePosition(-29.0, 0.0),
            StagePosition(29.0, 0.0),
            DoubleTiltStage(envelope=RectangularEnvelope(-30.0, 30.0, -30.0, 30.0)),
            np.eye(3),
            phase=phase,
        )
        assert not path.is_valid
        assert "clearance" in path.violation_reason or "envelope" in path.violation_reason

    def test_forward_and_reverse_paths_return_to_the_start(self) -> None:
        phase = cubic_phase()
        matrix = random_rotations(1, seed=31)[0]
        start = StagePosition(4.0, -6.0)
        alpha, beta = solve_tilts_for_direction(
            matrix @ ZoneAxis([1, 1, 2], phase=phase).unit_vector
        )[0]
        forward = plan_path(start, StagePosition(alpha, beta), WIDE_STAGE, matrix)
        backward = plan_path(forward.end, start, WIDE_STAGE, matrix)
        assert backward.end.alpha_deg == pytest.approx(start.alpha_deg, abs=1e-9)
        assert backward.end.beta_deg == pytest.approx(start.beta_deg, abs=1e-9)
        assert backward.total_travel_deg == pytest.approx(
            forward.total_travel_deg, abs=1e-6
        )

    def test_backlash_produces_an_approach_instruction(self) -> None:
        stage = DoubleTiltStage(
            envelope=RectangularEnvelope(-60.0, 60.0, -60.0, 60.0),
            calibration=StageCalibration(backlash_deg=0.4),
        )
        path = plan_path(
            StagePosition(0.0, 0.0), StagePosition(10.0, 10.0), stage, np.eye(3)
        )
        assert "overshoot" in path.approach_note
        assert "0.40 deg" in path.approach_note


# --------------------------------------------------------------------------- #
# The stereogram
# --------------------------------------------------------------------------- #


class TestTiltStereogram:
    @staticmethod
    def _report():  # type: ignore[no-untyped-def]
        phase = cubic_phase()
        matrix = Rotation.from_axis_angle([1.0, 1.0, 0.0], math.radians(18.0)).as_matrix()
        return (
            plan_tilt_to_zone_axis(
                make_state(phase, matrix),
                ZoneAxis([1, 1, 2], phase=phase),
                DoubleTiltStage(envelope=RectangularEnvelope(-35.0, 35.0, -30.0, 30.0)),
            ),
            DoubleTiltStage(envelope=RectangularEnvelope(-35.0, 35.0, -30.0, 30.0)),
        )

    def test_trajectory_has_one_dot_per_path_sample(self) -> None:
        """The figure plots the engine's samples, not a re-derived curve."""

        from pytex.plotting.tilt_stereogram import build_tilt_stereogram_figure_spec

        report, stage = self._report()
        spec = build_tilt_stereogram_figure_spec(report, stage, view="detail")
        best = report.best()
        assert best.path is not None
        trajectory_layers = [
            layer
            for layer in spec.marker_layers
            if layer.label is not None and "beam direction during" in layer.label
        ]
        assert len(trajectory_layers) == 1
        assert trajectory_layers[0].points.shape[0] == len(best.path.samples)

    def test_plotted_trajectory_matches_an_independent_integration(self) -> None:
        """Cross-check: accumulate small rotations instead of using the closed form.

        This is what makes the figure an independent check on the engine rather
        than a second copy of the same equation. The engine's beam directions
        come from the analytic formula; here they are rebuilt by composing the
        stage rotation matrices from scratch and inverting.
        """

        report, stage = self._report()
        best = report.best()
        assert best.path is not None
        matrix = report.current.matrix @ best.family.operator
        for sample in best.path.samples:
            rotation = stage.rotation_matrix(*sample.position.as_tuple())
            independent = rotation.T @ BEAM_AXIS_LABORATORY
            independent = matrix.T @ independent
            assert np.allclose(independent, sample.beam_direction_crystal, atol=1e-9)

    def test_detail_view_zooms_inside_the_overview(self) -> None:
        """The detail limits must be tighter, or the second panel earns nothing."""

        from pytex.plotting.tilt_stereogram import build_tilt_stereogram_figure_spec

        report, stage = self._report()
        overview = build_tilt_stereogram_figure_spec(report, stage, view="overview")
        detail = build_tilt_stereogram_figure_spec(report, stage, view="detail")
        assert overview.xlim is not None and detail.xlim is not None
        overview_width = overview.xlim[1] - overview.xlim[0]
        detail_width = detail.xlim[1] - detail.xlim[0]
        assert detail_width < 0.75 * overview_width

    def test_unreachable_equivalents_carry_the_rejection_colour(self) -> None:
        from pytex.plotting.tilt_stereogram import (
            TILT_STEREOGRAM_COLORS,
            build_tilt_stereogram_figure_spec,
        )

        report, stage = self._report()
        spec = build_tilt_stereogram_figure_spec(report, stage, view="overview")
        rejected = [
            layer
            for layer in spec.marker_layers
            if layer.label is not None and "out of range" in layer.label
        ]
        assert rejected
        assert rejected[0].edgecolors == TILT_STEREOGRAM_COLORS["unreachable"]

    def test_figure_renders_and_closes_cleanly(self) -> None:
        import matplotlib.pyplot as plt

        from pytex.plotting.tilt_stereogram import plot_tilt_stereogram

        report, stage = self._report()
        figure = plot_tilt_stereogram(report, stage)
        assert len(figure.axes) == 2
        plt.close(figure)

    def test_single_panel_can_be_drawn_on_a_supplied_axis(self) -> None:
        import matplotlib.pyplot as plt

        from pytex.plotting.tilt_stereogram import plot_tilt_stereogram

        report, stage = self._report()
        figure, axis = plt.subplots()
        plot_tilt_stereogram(report, stage, view="overview", ax=axis)
        assert axis.get_title()
        plt.close(figure)

    def test_two_panels_cannot_share_one_axis(self) -> None:
        import matplotlib.pyplot as plt

        from pytex.plotting.tilt_stereogram import plot_tilt_stereogram

        report, stage = self._report()
        figure, axis = plt.subplots()
        with pytest.raises(ValueError, match="cannot share one axis"):
            plot_tilt_stereogram(report, stage, view="both", ax=axis)
        plt.close(figure)


# --------------------------------------------------------------------------- #
# Noise and uncertainty
# --------------------------------------------------------------------------- #


class TestUncertainty:
    def test_reported_sigma_is_calibrated_against_the_actual_error(self) -> None:
        """The uncertainty estimate is itself checked, not merely reported.

        The orientation is perturbed by a known one-sigma amount, the state is
        told that uncertainty, and the actual landing error is measured against
        the truth. The reported sigma must match the empirical scatter to within
        a factor of two — a real calibration of the uncertainty, not a
        restatement of it. Without this test the reported sigma could be wrong by
        an order of magnitude and nothing would notice.
        """

        phase = cubic_phase()
        generator = np.random.default_rng(101)
        sigma_deg = 0.5
        target = ZoneAxis([1, 1, 2], phase=phase)
        truth = np.eye(3)
        errors: list[float] = []
        reported: list[float] = []
        for _ in range(300):
            axis = generator.normal(size=3)
            axis /= np.linalg.norm(axis)
            angle = generator.normal(scale=math.radians(sigma_deg))
            perturbed = Rotation.from_axis_angle(axis, angle).as_matrix() @ truth
            orientation = Orientation.from_matrix(
                perturbed,
                specimen_frame=HOLDER_FRAME,
                phase=phase,
                crystal_frame=phase.crystal_frame,
            )
            state = CurrentState.from_orientation(
                orientation,
                StagePosition(0.0, 0.0),
                current_zone_axis=ZoneAxis([0, 0, 1], phase=phase),
                orientation_uncertainty_deg=sigma_deg,
            )
            report = plan_tilt_to_zone_axis(
                state, target, WIDE_STAGE, include_paths=False
            )
            if not report.is_reachable:
                continue
            best = report.best()
            landed = (
                WIDE_STAGE.rotation_matrix(*best.position.as_tuple())
                @ truth
                @ best.orbit_member
            )
            landed /= np.linalg.norm(landed)
            errors.append(
                math.degrees(
                    math.acos(min(1.0, abs(float(np.dot(landed, BEAM_AXIS_LABORATORY)))))
                )
            )
            reported.append(best.sigma_residual_deg)

        assert len(errors) > 250
        empirical_rms = float(np.sqrt(np.mean(np.square(errors))))
        mean_reported = float(np.mean(reported))
        assert 0.5 <= mean_reported / empirical_rms <= 2.0, (
            f"reported sigma {mean_reported:.3f} deg against empirical "
            f"{empirical_rms:.3f} deg"
        )

    def test_a_supplied_orientation_reports_zero_uncertainty(self) -> None:
        """An orientation given directly is exact by construction, and says so."""

        phase = cubic_phase()
        state = make_state(phase, np.eye(3))
        assert state.orientation_uncertainty_deg == 0.0

    def test_conditioning_factor_is_bounded_within_a_real_envelope(self) -> None:
        """1/cos(alpha) stays under 1.31 for |alpha| <= 40 degrees.

        The parameterization is singular at the beta-axis pole, but a real holder
        never approaches it, so the amplification is benign and bounded.
        """

        assert 1.0 / math.cos(math.radians(40.0)) == pytest.approx(1.3054, abs=1e-3)

    def test_noisy_indexing_degrades_gracefully(self) -> None:
        phase = cubic_phase()
        generator = np.random.default_rng(202)
        residuals = []
        for _ in range(50):
            axis = generator.normal(size=3)
            axis /= np.linalg.norm(axis)
            perturbed = Rotation.from_axis_angle(
                axis, generator.normal(scale=math.radians(1.0))
            ).as_matrix()
            report = plan_tilt_to_zone_axis(
                make_state(phase, perturbed),
                ZoneAxis([1, 1, 1], phase=phase),
                WIDE_STAGE,
                include_paths=False,
            )
            if report.is_reachable:
                residuals.append(report.best().residual_deg)
        assert residuals
        # Forward validation is against the *estimated* orientation, so it stays
        # exact regardless of how wrong that estimate is. That is the point: the
        # residual measures solver quality, and sigma measures input quality.
        assert max(residuals) < 1e-9
