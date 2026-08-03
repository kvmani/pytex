"""TX5: solving a measured SAED pattern from picked spots.

Expected values are analytic identities of the geometry — the cubic ratio and
angle relations, the reciprocal relation d = 1/|g|, and closure of the
simulate-then-solve round trip — never a stored indexing.
"""

from __future__ import annotations

import json

import numpy as np
import pytest
from numpy.testing import assert_allclose

from pytex.adapters import measured_saed_pattern_schema_path
from pytex.core import (
    OrientationRelationship,
    Phase,
    Rotation,
    ZoneAxis,
)
from pytex.diffraction import (
    MEASURED_PATTERN_SCHEMA,
    PATTERN_SOLUTION_SCHEMA,
    KinematicSimulationConfig,
    MeasuredSAEDPattern,
    MeasuredSpot,
    PatternCalibration,
    SpotTable,
    assign_transformation_variant,
    simulate_composite_saed,
    simulate_composite_saed_from_child_zone,
    simulate_zone_axis_spots,
    solve_saed_pattern,
    solve_saed_pattern_file,
)
from pytex.plotting.saed_picker import SpotPickerState
from tests.unit.test_composite_saed import make_bcc_hcp_phases, make_fcc_bcc_phases

CAMERA_CONSTANT_MM_ANGSTROM = 180.0


def _config() -> KinematicSimulationConfig:
    return KinematicSimulationConfig(
        camera_constant_mm_angstrom=CAMERA_CONSTANT_MM_ANGSTROM
    )


def _pattern_from_simulation(
    phase: Phase,
    zone_indices: tuple[int, int, int],
    *,
    name: str = "synthetic",
    noise_mm: float = 0.0,
    seed: int = 0,
) -> tuple[MeasuredSAEDPattern, SpotTable]:
    """Turn a simulated pattern into a measured one, as a user's picks would be."""

    zone = ZoneAxis(np.asarray(zone_indices, dtype=np.int64), phase=phase)
    spots = simulate_zone_axis_spots(phase, zone, config=_config())
    coordinates = np.asarray(spots.detector_mm, dtype=np.float64)
    if noise_mm > 0.0:
        rng = np.random.default_rng(seed)
        coordinates = coordinates + rng.normal(0.0, noise_mm, size=coordinates.shape)
    measured = MeasuredSAEDPattern(
        name=name,
        spots=tuple(
            MeasuredSpot(position=(float(x), float(y))) for x, y in coordinates
        ),
        calibration=PatternCalibration(
            units="mm", camera_constant_mm_angstrom=CAMERA_CONSTANT_MM_ANGSTROM
        ),
    )
    return measured, spots


class TestCalibration:
    def test_millimetres_convert_through_the_camera_constant(self) -> None:
        """|g| = r / (L*lambda) is the definition of the camera constant."""

        calibration = PatternCalibration(units="mm", camera_constant_mm_angstrom=200.0)
        converted = calibration.to_reciprocal_angstrom([[100.0, 0.0], [0.0, -50.0]])
        assert_allclose(converted, [[0.5, 0.0], [0.0, -0.25]], rtol=1e-12)

    def test_pixels_scale_through_the_pixel_size_then_the_camera_constant(self) -> None:
        calibration = PatternCalibration(
            units="px",
            pixel_size_mm=0.02,
            camera_constant_mm_angstrom=200.0,
            centre=(512.0, 512.0),
        )
        # 100 px from the centre is 2 mm, hence 0.01 1/angstrom.
        converted = calibration.to_reciprocal_angstrom([[612.0, 512.0]])
        assert_allclose(converted, [[0.01, 0.0]], rtol=1e-12)

    def test_camera_length_and_voltage_reproduce_the_camera_constant(self) -> None:
        """L*lambda with the relativistic 200 kV wavelength, 0.02508 angstrom."""

        calibration = PatternCalibration(
            units="mm", camera_length_mm=800.0, beam_energy_kev=200.0
        )
        assert_allclose(
            calibration.effective_camera_constant_mm_angstrom, 800.0 * 0.02508, rtol=2e-4
        )

    def test_reciprocal_units_pass_through_unchanged(self) -> None:
        calibration = PatternCalibration(units="reciprocal_angstrom")
        assert calibration.effective_camera_constant_mm_angstrom is None
        assert_allclose(
            calibration.to_reciprocal_angstrom([[0.3, -0.4]]), [[0.3, -0.4]], atol=0.0
        )

    def test_pixels_without_a_pixel_size_are_rejected_at_construction(self) -> None:
        with pytest.raises(ValueError, match="positive pixel_size_mm"):
            PatternCalibration(units="px", camera_constant_mm_angstrom=200.0)

    def test_lengths_without_a_camera_constant_are_rejected_at_construction(self) -> None:
        """An uncalibrated length is unrecoverable, so fail early, not per spot."""

        with pytest.raises(ValueError, match="need a camera constant"):
            PatternCalibration(units="mm")

    def test_unknown_units_are_rejected(self) -> None:
        with pytest.raises(ValueError, match="units must be one of"):
            PatternCalibration(units="nanometres")


class TestMeasuredPattern:
    def test_observed_d_spacings_are_the_reciprocal_of_g(self) -> None:
        phase, _ = make_fcc_bcc_phases()
        measured, _ = _pattern_from_simulation(phase, (0, 0, 1))
        assert_allclose(
            measured.d_spacings_angstrom(),
            1.0 / measured.g_magnitudes_inv_angstrom(),
            rtol=1e-12,
        )

    def test_a_pattern_needs_two_spots(self) -> None:
        with pytest.raises(ValueError, match="at least two spots"):
            MeasuredSAEDPattern(
                name="one",
                spots=(MeasuredSpot(position=(1.0, 0.0)),),
                calibration=PatternCalibration(),
            )

    def test_a_spot_on_the_transmitted_beam_is_rejected(self) -> None:
        """The beam is the calibration centre, not a reflection."""

        with pytest.raises(ValueError, match="coincides with the transmitted beam"):
            MeasuredSAEDPattern(
                name="beam",
                spots=(
                    MeasuredSpot(position=(0.0, 0.0)),
                    MeasuredSpot(position=(1.0, 0.0)),
                ),
                calibration=PatternCalibration(),
            )

    def test_yaml_round_trip_preserves_every_coordinate(self, tmp_path) -> None:
        phase, _ = make_fcc_bcc_phases()
        measured, _ = _pattern_from_simulation(phase, (1, 1, 0), name="round_trip")
        path = measured.to_yaml(tmp_path / "pattern.yaml")
        restored = MeasuredSAEDPattern.from_yaml(path)
        assert restored.name == measured.name
        assert len(restored) == len(measured)
        assert_allclose(
            restored.g_vectors_inv_angstrom(),
            measured.g_vectors_inv_angstrom(),
            atol=0.0,
        )

    def test_yaml_payload_validates_against_the_schema(self, tmp_path) -> None:
        jsonschema = pytest.importorskip("jsonschema")
        phase, _ = make_fcc_bcc_phases()
        measured, _ = _pattern_from_simulation(phase, (1, 1, 1))
        payload = measured.to_dict()
        assert payload["schema"] == MEASURED_PATTERN_SCHEMA
        jsonschema.validate(
            json.loads(json.dumps(payload)),
            json.loads(measured_saed_pattern_schema_path().read_text(encoding="utf-8")),
        )

    def test_an_unknown_schema_is_refused(self) -> None:
        with pytest.raises(ValueError, match="Unsupported measured-pattern schema"):
            MeasuredSAEDPattern.from_dict(
                {
                    "schema": "something.else/9",
                    "name": "x",
                    "calibration": {"units": "reciprocal_angstrom"},
                    "spots": [{"x": 1.0, "y": 0.0}, {"x": 0.0, "y": 1.0}],
                }
            )


class TestSyntheticClosure:
    """Simulate, then solve: the answer must be the pattern that was simulated."""

    @pytest.mark.parametrize(
        "zone", [(0, 0, 1), (0, 1, 1), (1, 1, 1), (1, 1, 2), (0, 1, 3)]
    )
    def test_every_spot_is_indexed_and_the_zone_axis_is_recovered(
        self, zone: tuple[int, int, int]
    ) -> None:
        parent, child = make_fcc_bcc_phases()
        measured, _ = _pattern_from_simulation(parent, zone)
        report = solve_saed_pattern(measured, [parent, child], max_index=6)
        best = report.best()
        assert best.phase_name == parent.name
        assert best.matched_fraction == 1.0
        assert best.mean_residual_inv_angstrom < 1e-9
        # The recovered zone axis must be the simulated one up to sign and to
        # the crystal symmetry, which is what the canonical description reports.
        recovered = np.abs(np.asarray(best.zone_axis.indices, dtype=np.int64))
        assert sorted(recovered.tolist()) == sorted(np.abs(np.asarray(zone)).tolist())

    def test_the_solved_indices_agree_with_the_simulated_ones_up_to_symmetry(self) -> None:
        """Index *families* must match spot by spot, whatever description is chosen."""

        parent, _ = make_fcc_bcc_phases()
        measured, spots = _pattern_from_simulation(parent, (0, 0, 1))
        report = solve_saed_pattern(measured, [parent], max_index=6)
        best = report.best()
        for solved in best.solved_spots:
            simulated = np.abs(spots.hkl[solved.measured_index])
            assert sorted(np.abs(solved.hkl).tolist()) == sorted(simulated.tolist())

    def test_a_hexagonal_phase_is_solved_and_labeled_in_four_index_form(self) -> None:
        _, alpha = make_bcc_hcp_phases()
        measured, _ = _pattern_from_simulation(alpha, (0, 0, 1))
        report = solve_saed_pattern(measured, [alpha], max_index=6)
        best = report.best()
        assert best.matched_fraction == 1.0
        assert all(spot.label.count(" ") == 3 for spot in best.solved_spots)


class TestCubicIdentities:
    def test_the_cube_zone_shows_the_root_two_ratio_and_forty_five_degree_angle(
        self,
    ) -> None:
        """Down [001] of an fcc crystal, {220} sits at sqrt(2) times {200}, at 45 deg.

        Independent geometry, not program output: |g_220| / |g_200| = sqrt(2)
        exactly, and the angle between (200) and (220) is exactly 45 degrees.
        The solver's own output must reproduce both.
        """

        parent, _ = make_fcc_bcc_phases()
        measured, _ = _pattern_from_simulation(parent, (0, 0, 1))
        report = solve_saed_pattern(measured, [parent], max_index=6)
        best = report.best()
        by_family: dict[tuple[int, ...], np.ndarray] = {}
        for spot in best.solved_spots:
            key = tuple(sorted(np.abs(spot.hkl).tolist()))
            by_family.setdefault(key, np.asarray(spot.predicted_g_inv_angstrom))
        two_hundred = by_family[(0, 0, 2)]
        two_twenty = by_family[(0, 2, 2)]
        ratio = float(np.linalg.norm(two_twenty) / np.linalg.norm(two_hundred))
        assert_allclose(ratio, np.sqrt(2.0), rtol=1e-9)
        cosine = float(
            np.dot(two_hundred, two_twenty)
            / (np.linalg.norm(two_hundred) * np.linalg.norm(two_twenty))
        )
        assert_allclose(abs(cosine), np.cos(np.radians(45.0)), atol=1e-9)


class TestDiscrimination:
    def test_an_fcc_pattern_is_not_solved_as_bcc(self) -> None:
        """Systematic absences differ, so the wrong lattice cannot explain the ratios."""

        parent, child = make_fcc_bcc_phases()
        measured, _ = _pattern_from_simulation(parent, (0, 0, 1))
        report = solve_saed_pattern(measured, [child], max_index=6)
        if report.solutions:
            assert report.best().matched_fraction < 1.0
        assert not report.is_conclusive

    def test_the_true_phase_outranks_the_decoy_when_both_are_offered(self) -> None:
        parent, child = make_fcc_bcc_phases()
        measured, _ = _pattern_from_simulation(parent, (0, 1, 1))
        report = solve_saed_pattern(measured, [child, parent], max_index=6)
        assert report.best().phase_name == parent.name
        assert report.is_conclusive

    def test_an_unsolvable_pattern_reports_no_solution_rather_than_guessing(self) -> None:
        """Two arbitrary spots that no candidate phase can explain."""

        parent, _ = make_fcc_bcc_phases()
        measured = MeasuredSAEDPattern(
            name="nonsense",
            spots=(
                MeasuredSpot(position=(0.001, 0.0)),
                MeasuredSpot(position=(0.0, 0.0013)),
            ),
            calibration=PatternCalibration(units="reciprocal_angstrom"),
        )
        report = solve_saed_pattern(measured, [parent], max_index=4)
        assert len(report) == 0
        assert not report.is_conclusive
        assert "could not be solved" in report.describe()
        with pytest.raises(ValueError, match="No solution was found"):
            report.best()


class TestNoiseRobustness:
    def test_moderate_picking_noise_still_recovers_the_pattern(self) -> None:
        parent, child = make_fcc_bcc_phases()
        measured, _ = _pattern_from_simulation(
            parent, (0, 1, 1), noise_mm=0.5, seed=3
        )
        report = solve_saed_pattern(measured, [parent, child], max_index=6)
        best = report.best()
        assert best.phase_name == parent.name
        assert best.matched_fraction == 1.0
        # The residual must reflect the noise rather than being absorbed silently.
        assert best.mean_residual_inv_angstrom > 0.0

    def test_severe_noise_degrades_to_a_partial_or_absent_solution(self) -> None:
        """The failure mode must be visible, not a confident wrong answer."""

        parent, child = make_fcc_bcc_phases()
        measured, _ = _pattern_from_simulation(
            parent, (0, 1, 1), noise_mm=25.0, seed=5
        )
        report = solve_saed_pattern(measured, [parent, child], max_index=6)
        assert not report.is_conclusive or report.best().matched_fraction < 1.0


class TestReportContract:
    def test_describe_names_the_intrinsic_zone_sense_ambiguity(self) -> None:
        parent, _ = make_fcc_bcc_phases()
        measured, _ = _pattern_from_simulation(parent, (1, 1, 1))
        report = solve_saed_pattern(measured, [parent], max_index=6)
        text = report.describe()
        assert "cannot distinguish a zone axis from its reverse" in text

    def test_json_payload_declares_its_schema_and_matches_the_solutions(self) -> None:
        parent, _ = make_fcc_bcc_phases()
        measured, _ = _pattern_from_simulation(parent, (0, 0, 1))
        report = solve_saed_pattern(measured, [parent], max_index=6)
        payload = json.loads(json.dumps(report.to_json_dict()))
        assert payload["schema"] == PATTERN_SOLUTION_SCHEMA
        assert payload["is_conclusive"] is report.is_conclusive
        assert len(payload["solutions"]) == len(report)
        assert len(payload["solutions"][0]["spots"]) == len(
            report.best().solved_spots
        )

    def test_unindexed_spots_are_reported_with_the_max_index_hint(self) -> None:
        """A spot beyond the solver's index bound is never offered a match."""

        parent, _ = make_fcc_bcc_phases()
        measured, _ = _pattern_from_simulation(parent, (0, 1, 1))
        report = solve_saed_pattern(measured, [parent], max_index=4)
        best = report.best()
        assert best.unindexed_spot_indices
        assert "max_index" in best.describe()

    def test_solving_from_a_file_matches_solving_in_memory(self, tmp_path) -> None:
        parent, _ = make_fcc_bcc_phases()
        measured, _ = _pattern_from_simulation(parent, (0, 0, 1), name="from_file")
        path = measured.to_yaml(tmp_path / "pattern.yaml")
        in_memory = solve_saed_pattern(measured, [parent], max_index=6)
        from_file = solve_saed_pattern_file(path, [parent], max_index=6)
        assert from_file.best().zone_axis_label == in_memory.best().zone_axis_label
        assert from_file.best().matched_fraction == in_memory.best().matched_fraction


class TestVariantAssignment:
    def test_a_solved_child_pattern_is_assigned_its_planted_variant(self) -> None:
        """The composite engine plants the variant; the solver must name it back.

        The pattern is anchored on the chosen variant's own basal zone, so its
        spots really do form a zone-axis pattern — which is what the solver
        assumes, and what an operator produces by tilting the product on zone.
        """

        beta, alpha = make_bcc_hcp_phases()
        relationship = OrientationRelationship.from_burgers_correspondence(
            parent_phase=beta, child_phase=alpha
        )
        planted = 3
        composite = simulate_composite_saed_from_child_zone(
            relationship,
            ZoneAxis(np.array([0, 0, 1]), phase=alpha),
            anchor_variant_index=planted,
            variant_indices=(planted,),
            config=_config(),
        )
        variant_pattern = composite.variant_pattern(planted)
        measured = MeasuredSAEDPattern(
            name="child_variant",
            spots=tuple(
                MeasuredSpot(position=(float(x), float(y)))
                for x, y in variant_pattern.spots.detector_mm
            ),
            calibration=PatternCalibration(
                units="mm", camera_constant_mm_angstrom=CAMERA_CONSTANT_MM_ANGSTROM
            ),
        )
        report = solve_saed_pattern(measured, [alpha], max_index=6)
        best = report.best()
        assert best.matched_fraction == 1.0
        # The shared detector basis has columns (u, v, z) in *parent crystal*
        # Cartesian coordinates, so a parent vector's pattern coordinates are
        # basis^T p: the parent's crystal-to-pattern rotation is basis^T.
        parent_orientation = Rotation.from_matrix(
            np.ascontiguousarray(composite.zone_basis_parent.T)
        )
        assigned = assign_transformation_variant(best, relationship, parent_orientation)
        assert assigned.variant_index == planted
        assert assigned.variant_deviation_deg is not None
        assert assigned.variant_deviation_deg < 1.0
        assert "variant" in assigned.describe()

    def test_an_off_zone_variant_pattern_is_only_partly_indexed(self) -> None:
        """A stated limitation, pinned rather than left to be discovered.

        The solver assumes the spots form a zone-axis pattern about a low-index
        zone. A variant seen from a *parent* zone axis is generally off its own
        zone — its child zone axis is irrational — so its excitation-selected
        spots do not all lie in one zero-order Laue zone and cannot all be
        indexed. The partial match is the honest outcome; a full match here
        would mean the solver was inventing reflections.
        """

        beta, alpha = make_bcc_hcp_phases()
        relationship = OrientationRelationship.from_burgers_correspondence(
            parent_phase=beta, child_phase=alpha
        )
        composite = simulate_composite_saed(
            relationship,
            ZoneAxis(np.array([1, 1, 0]), phase=beta),
            variant_indices=(3,),
            config=_config(),
        )
        variant_pattern = composite.variant_pattern(3)
        assert variant_pattern.nearest_zone_axis.deviation_deg > 1.0
        measured = MeasuredSAEDPattern(
            name="off_zone_variant",
            spots=tuple(
                MeasuredSpot(position=(float(x), float(y)))
                for x, y in variant_pattern.spots.detector_mm
            ),
            calibration=PatternCalibration(
                units="mm", camera_constant_mm_angstrom=CAMERA_CONSTANT_MM_ANGSTROM
            ),
        )
        report = solve_saed_pattern(measured, [alpha], max_index=6)
        assert not report.is_conclusive or report.best().matched_fraction < 1.0

    def test_a_mismatched_child_phase_is_rejected(self) -> None:
        beta, alpha = make_bcc_hcp_phases()
        relationship = OrientationRelationship.from_burgers_correspondence(
            parent_phase=beta, child_phase=alpha
        )
        measured, _ = _pattern_from_simulation(beta, (0, 0, 1))
        report = solve_saed_pattern(measured, [beta], max_index=6)
        with pytest.raises(ValueError, match="child phase"):
            assign_transformation_variant(
                report.best(), relationship, Rotation.identity()
            )


class TestPickerStateMachine:
    """The picking logic is a plain object, so it is tested without a display."""

    def test_add_remove_undo_and_clear(self) -> None:
        state = SpotPickerState()
        state.add(1.0, 0.0)
        state.add(0.0, 1.0)
        state.add(2.0, 2.0)
        assert len(state) == 3
        assert state.remove_nearest(0.1, 0.1) == 0
        assert len(state) == 2
        assert state.undo() == 1
        assert len(state) == 1
        state.clear()
        assert len(state) == 0

    def test_removal_radius_prevents_deleting_a_distant_pick(self) -> None:
        state = SpotPickerState()
        state.add(10.0, 10.0)
        assert state.remove_nearest(0.0, 0.0, radius=1.0) is None
        assert len(state) == 1

    def test_undo_and_remove_on_an_empty_session_do_nothing(self) -> None:
        state = SpotPickerState()
        assert state.undo() is None
        assert state.remove_nearest(0.0, 0.0) is None

    def test_the_session_centre_overrides_the_calibrations(self) -> None:
        """The centre is something the user clicks, so their click wins."""

        state = SpotPickerState()
        state.add(3.0, 0.0)
        state.add(0.0, 4.0)
        state.set_centre(1.0, 1.0)
        pattern = state.to_pattern(
            name="picked",
            calibration=PatternCalibration(
                units="mm", camera_constant_mm_angstrom=100.0, centre=(99.0, 99.0)
            ),
        )
        assert pattern.calibration.centre == (1.0, 1.0)
        assert_allclose(
            pattern.g_vectors_inv_angstrom(), [[0.02, -0.01], [-0.01, 0.03]], rtol=1e-12
        )

    def test_a_picked_session_round_trips_through_yaml_and_solves(self, tmp_path) -> None:
        parent, _ = make_fcc_bcc_phases()
        _, spots = _pattern_from_simulation(parent, (0, 0, 1))
        state = SpotPickerState()
        for x, y in spots.detector_mm:
            state.add(float(x), float(y))
        pattern = state.to_pattern(
            name="picked_cube_zone",
            calibration=PatternCalibration(
                units="mm", camera_constant_mm_angstrom=CAMERA_CONSTANT_MM_ANGSTROM
            ),
        )
        path = pattern.to_yaml(tmp_path / "picked.yaml")
        report = solve_saed_pattern_file(path, [parent], max_index=6)
        assert report.best().matched_fraction == 1.0

    def test_non_finite_picks_are_rejected(self) -> None:
        state = SpotPickerState()
        with pytest.raises(ValueError, match="must be finite"):
            state.add(float("nan"), 0.0)
        with pytest.raises(ValueError, match="must be finite"):
            state.set_centre(0.0, float("inf"))


class TestValidation:
    def test_no_candidate_phases_is_rejected(self) -> None:
        parent, _ = make_fcc_bcc_phases()
        measured, _ = _pattern_from_simulation(parent, (0, 0, 1))
        with pytest.raises(ValueError, match="at least one candidate phase"):
            solve_saed_pattern(measured, [])

    def test_out_of_range_tolerances_are_rejected(self) -> None:
        parent, _ = make_fcc_bcc_phases()
        measured, _ = _pattern_from_simulation(parent, (0, 0, 1))
        with pytest.raises(ValueError, match="length_tolerance_relative"):
            solve_saed_pattern(measured, [parent], length_tolerance_relative=1.5)
        with pytest.raises(ValueError, match="angle_tolerance_deg"):
            solve_saed_pattern(measured, [parent], angle_tolerance_deg=0.0)
